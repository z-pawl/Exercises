classdef FlowComponent < IComponent
    properties
        noi (1,1)           % Number of iterations

        % Dependencies
        grid Grid2D {mustBeScalarOrEmpty}


        % Collocated grid approach is used to represent velocity and pressure fields 
        vx (:,:) double             % X - component of velocity
        vr (:,:) double             % R - component of velocity
        p (:,:) double              % Pressure

        % Rhie-chow interpolated velocities
        vx_faces (:,:) double
        vr_faces (:,:) double


        % Physical properties
        porosity double             % Porosity
        rho (:,:) double            % Density
        visc (:,:) double           % Dynamic viscosity

        % Source terms
        srcx (:,:) double           % Constant momentum source term in the x direction
        src_linx (:,:) double       % Linear momentum source term in the x direction
        srcr (:,:) double           % Constant momentum source term in the r direction
        src_linr (:,:) double       % Linear momentum source term in the r direction


        % Functions used to calculate physical properties
        rho_function (1,1) function_handle = @(x) ones(x.grid.sz);
        visc_function (1,1) function_handle = @(x) ones(x.grid.sz);
        
        % Functions used to calculate source terms
        srcx_function (1,1) function_handle = @(x) zeros(x.grid.sz);
        src_linx_function (1,1) function_handle = @(x) zeros(x.grid.sz);
        srcr_function (1,1) function_handle = @(x) zeros(x.grid.sz);
        src_linr_function (1,1) function_handle = @(x) zeros(x.grid.sz);


        % Boundary conditions
        vx_bds Boundaries {mustBeScalarOrEmpty}         % X - component velocity boundary conditions
        vr_bds Boundaries {mustBeScalarOrEmpty}         % R - component velocity boundary conditions
        p_bds Boundaries {mustBeScalarOrEmpty}          % Pressure boundary conditions
        p_corr_bds Boundaries {mustBeScalarOrEmpty}     % Pressure correction boundary conditions


        % Solver settings
        tol_v (1,1) double {mustBeNonnegative}          % Convergence criteria for velocity
        tol_p (1,1) double {mustBeNonnegative}          % Convergence criteria for pressure
        relaxation_factor_v (1,1) double                % Velocity relaxation factor
        relaxation_factor_p (1,1) double                % Pressure relaxation factor
        inner_iters (1,1)                               % Number of inner iterations per outer iteration
        solver_iters_v (1,1)                            % Number of velocity solver iterations per inner iteration
        solver_iters_p (1,1)                            % Number of pressure solver iterations per inner iteration
        residual_history_vx (:,1) double                % X component of velocity residual history
        residual_history_vr (:,1) double                % R component of velocity residual history
        residual_history_p (:,1) double                 % Pressure residual history
    end
    methods
        function obj = FlowComponent(grid, vx, vr, p, porosity, rho_function, visc_function, srcx_function, src_linx_function, srcr_function, src_linr_function, tol_v, tol_p, relaxation_factor_v, relaxation_factor_p, inner_iters, solver_iters_v, solver_iters_p)
            arguments
                grid (1,1) Grid2D
                vx (:,:) double     % Initial x-component velocity
                vr (:,:) double     % Initial r-component velocity
                p  (:,:) double     % Initial pressure
                porosity double
                rho_function (1,1) function_handle
                visc_function (1,1) function_handle
                srcx_function (1,1) function_handle
                src_linx_function (1,1) function_handle
                srcr_function (1,1) function_handle
                src_linr_function (1,1) function_handle
                tol_v (1,1) double
                tol_p (1,1) double
                relaxation_factor_v (1,1) double 
                relaxation_factor_p (1,1) double
                inner_iters (1,1)
                solver_iters_v (1,1)
                solver_iters_p (1,1)
            end

            if ~isequal(size(vx), grid.sz) || ~isequal(size(vr), grid.sz) || ~isequal(size(p), grid.sz)
                error("The initial fields must have the same size as the grid");
            end
            
            obj.noi = 0;

            obj.grid = grid;

            obj.vx = vx;
            obj.vr = vr;
            obj.p = p;

            obj.vx_faces = zeros(grid.sz + [1 0]);
            obj.vr_faces = zeros(grid.sz + [0 1]);

            obj.porosity = porosity;
            obj.rho_function = rho_function;
            obj.visc_function = visc_function;

            obj.srcx_function = srcx_function;
            obj.src_linx_function = src_linx_function;
            obj.srcr_function = srcr_function;
            obj.src_linr_function = src_linr_function;

            obj.vx_bds = Boundaries(grid);
            obj.vr_bds = Boundaries(grid);
            obj.p_bds = Boundaries(grid);
            obj.p_corr_bds = Boundaries(grid);

            obj.tol_v = tol_v;
            obj.tol_p = tol_p;
            obj.relaxation_factor_v = relaxation_factor_v;
            obj.relaxation_factor_p = relaxation_factor_p;
            obj.inner_iters = inner_iters;
            obj.solver_iters_v = solver_iters_v;
            obj.solver_iters_p = solver_iters_p;
            obj.residual_history_vx = zeros(10000,1);
            obj.residual_history_vr = zeros(10000,1);
            obj.residual_history_p = zeros(10000,1);
        end

        function update_properties(obj)
            obj.rho = obj.rho_function();
            obj.visc = obj.visc_function();

            obj.srcx = obj.srcx_function();
            obj.src_linx = obj.src_linx_function();
            obj.srcr = obj.srcr_function();
            obj.src_linr = obj.src_linr_function();
        end

        % Updates face velocities, face velocities are calculated using Rhie-Chow interpolation
        function update_face_velocities(obj, vx, vr, aP_vx, aP_vr)
            [rhie_chow_vx, rhie_chow_vr] = rhie_chow_interpolation(obj.grid, vx, vr, obj.p, obj.vx_bds, obj.vr_bds, obj.p_bds, aP_vx, aP_vr);
            
            obj.vx_faces = rhie_chow_vx;
            obj.vr_faces = rhie_chow_vr;
        end

        % Converts pressure boundary conditions to pressure correction boundary conditions
        function convert_pressure_bd_conditions(obj)
            obj.p_corr_bds.boundary_conditions_x.is_bd = obj.p_bds.boundary_conditions_x.is_bd;
            obj.p_corr_bds.boundary_conditions_x.orientation = obj.p_bds.boundary_conditions_x.orientation;
            obj.p_corr_bds.boundary_conditions_x.a = obj.p_bds.boundary_conditions_x.a;
            obj.p_corr_bds.boundary_conditions_x.b = obj.p_bds.boundary_conditions_x.b;
            obj.p_corr_bds.boundary_conditions_x.c = zeros(obj.grid.sz + [1 0]);

            obj.p_corr_bds.boundary_conditions_r.is_bd = obj.p_bds.boundary_conditions_r.is_bd;
            obj.p_corr_bds.boundary_conditions_r.orientation = obj.p_bds.boundary_conditions_r.orientation;
            obj.p_corr_bds.boundary_conditions_r.a = obj.p_bds.boundary_conditions_r.a;
            obj.p_corr_bds.boundary_conditions_r.b = obj.p_bds.boundary_conditions_r.b;
            obj.p_corr_bds.boundary_conditions_r.c = zeros(obj.grid.sz + [0 1]);
        end

        function [coeff_vx, coeff_vr] = get_coefficients_v(obj)
            % Coefficients used to calculate the physical properties at faces
            [lin_int_coeffs_x, lin_int_coeffs_r] = linear_interpolation_scheme(obj.grid, obj.grid.domain_boundary);
            
            % Properties interpolated to faces
            % Density
            [rho_x, rho_r] = evaluate_faces(lin_int_coeffs_x, lin_int_coeffs_r, obj.rho);
            % Dynamic viscosity
            [visc_x, visc_r] = evaluate_faces(lin_int_coeffs_x, lin_int_coeffs_r, obj.visc);

            clear lin_int_coeffs_x lin_int_coeffs_r

            % F = rho * v * A / porosity ^ 2
            Fx = rho_x .* obj.vx_faces .* obj.grid.face_area_x / obj.porosity .^ 2;
            Fr = rho_r .* obj.vr_faces .* obj.grid.face_area_r / obj.porosity .^ 2;

            % D = visc * A / porosity
            Dx = visc_x .* obj.grid.face_area_x / obj.porosity;
            Dr = visc_r .* obj.grid.face_area_r / obj.porosity;

            clear rho_x rho_r visc_x visc_r;


            % Calculating velocities coefficients
            coeff_vx = obj.get_coefficients_vx(Fx, Fr, Dx, Dr);
            coeff_vr = obj.get_coefficients_vr(Fx, Fr, Dx, Dr);
            clear Fx Fr Dx Dr;


            % Calculating the volume integral of pressure gradient using the Green Gauss method

            % Coefficients for pressure interpolation
            [coeff_px, coeff_pr] = linear_interpolation_scheme(obj.grid, obj.p_bds);
            % Pressure values at faces
            [px, pr] = evaluate_faces(coeff_px, coeff_pr, obj.p);
            % Multiplying the pressure by face areas
            px = px .* obj.grid.face_area_x;
            pr = pr .* obj.grid.face_area_r;
            % Adding the volume integral of pressure gradient as source terms
            coeff_vx(:,:,6) = coeff_vx(:,:,6) + px(1:end-1,:) - px(2:end,:);
            coeff_vr(:,:,6) = coeff_vr(:,:,6) + pr(:,1:end-1) - pr(:,2:end);
        end

        function coeff_vx = get_coefficients_vx(obj, Fx, Fr, Dx, Dr)
            % Values and normal derivatives of velocity at faces

            % Coefficients for normal derivatives
            [coeffs_x_n_der_vx, coeffs_r_n_der_vx] = central_differencing_scheme(obj.grid, obj.vx_bds);

            % Values of normal derivatives
            [n_der_vx_x, n_der_vx_r] = evaluate_faces(coeffs_x_n_der_vx, coeffs_r_n_der_vx, obj.vx);

            % Coefficients for the upwind scheme
            [coeffs_x_upwind, coeffs_r_upwind] = upwind_scheme(obj.grid, obj.vx_faces, obj.vr_faces, obj.vx_bds);
            
            % Coefficients for the TVD scheme
            [coeffs_x_TVD, coeffs_r_TVD] = TVD_scheme(obj.grid, obj.vx_faces, obj.vr_faces, @van_leer_flux_limiter, n_der_vx_x, n_der_vx_r, obj.vx_bds);

            % Coefficeints for the deferred scheme
            [coeffs_x_vx_deferred, coeffs_r_vx_deferred] = deferred_correction_face(coeffs_x_upwind, coeffs_r_upwind, coeffs_x_TVD, coeffs_r_TVD, obj.vx);

            clear n_der_vx_x n_der_vx_r coeffs_x_upwind coeffs_r_upwind coeffs_x_TVD coeffs_r_TVD;

            coeff_vx = assemble_coeff_array(Fx, Fr, Dx, Dr, coeffs_x_vx_deferred, coeffs_r_vx_deferred, coeffs_x_n_der_vx, coeffs_r_n_der_vx, obj.srcx, obj.src_linx, obj.grid.volume);
        end

        function coeff_vr = get_coefficients_vr(obj, Fx, Fr, Dx, Dr)
            % Values and normal derivatives of velocity at faces

            % Coefficients for normal derivatives
            [coeffs_x_n_der_vr, coeffs_r_n_der_vr] = central_differencing_scheme(obj.grid, obj.vr_bds);

            % Values of normal derivatives
            [n_der_vr_x, n_der_vr_r] = evaluate_faces(coeffs_x_n_der_vr, coeffs_r_n_der_vr, obj.vr);

            % Coefficients for the upwind scheme
            [coeffs_x_upwind, coeffs_r_upwind] = upwind_scheme(obj.grid, obj.vx_faces, obj.vr_faces, obj.vr_bds);
            
            % Coefficients for the TVD scheme
            [coeffs_x_TVD, coeffs_r_TVD] = TVD_scheme(obj.grid, obj.vx_faces, obj.vr_faces, @van_leer_flux_limiter, n_der_vr_x, n_der_vr_r, obj.vr_bds);

            % Coefficeints for the deferred scheme
            [coeffs_x_vr_deferred, coeffs_r_vr_deferred] = deferred_correction_face(coeffs_x_upwind, coeffs_r_upwind, coeffs_x_TVD, coeffs_r_TVD, obj.vr);

            clear n_der_vr_x n_der_vr_r coeffs_x_upwind coeffs_r_upwind coeffs_x_TVD coeffs_r_TVD;

            coeff_vr = assemble_coeff_array(Fx, Fr, Dx, Dr, coeffs_x_vr_deferred, coeffs_r_vr_deferred, coeffs_x_n_der_vr, coeffs_r_n_der_vr, obj.srcr, obj.src_linr, obj.grid.volume);
        
            % Additional geometric term
            coeff_vr(:,:,1) = coeff_vr(:,:,1) + obj.visc / (obj.porosity .* obj.grid.cent_pos_r .^ 2);
        end

        function coeff_p_corr = get_coefficients_p_corr(obj, coeff_vx_aP, coeff_vr_aP)
            % Coefficients used to calculate the linear interpolation
            [lin_int_coeffs_x, lin_int_coeffs_r] = linear_interpolation_scheme(obj.grid, obj.grid.domain_boundary);

            % Interpolated quantities
            % Density
            [rho_x, rho_r] = evaluate_faces(lin_int_coeffs_x, lin_int_coeffs_r, obj.rho);
            % Central momentum coefficients
            [coeff_vx_aP_f, ~] = evaluate_faces(lin_int_coeffs_x, lin_int_coeffs_r, coeff_vx_aP);
            [~, coeff_vr_aP_f] = evaluate_faces(lin_int_coeffs_x, lin_int_coeffs_r, coeff_vr_aP);

            clear lin_int_coeffs_x lin_int_coeffs_x coeff_vx_aP coeff_vr_aP

            % Calculating the coefficients used in the discretized equations
            c_x = rho_x .* obj.grid.face_area_x ./ coeff_vx_aP_f;
            c_r = rho_r .* obj.grid.face_area_r ./ coeff_vr_aP_f;
            Fx = rho_x .* obj.vx_faces .* obj.grid.face_area_x;
            Fr = rho_r .* obj.vr_faces .* obj.grid.face_area_r;
            clear rho_x rho_r coeff_vx_aP_f coeff_vr_aP_f

            %%%% DISCRETIZATION SCHEME FOR CALCULATING THE PRESSURE DERIVATIVE AT A FACE HAS TO BE CHANGED TO A GREEN GAUSS BASED METHOD
            % Calculating coefficients for normal derivatives at faces for pressure
            [coeff_der_p_x, coeff_der_p_r] = central_differencing_scheme(obj.grid, obj.p_corr_bds);


            % Preallocating the array
            coeff_p_corr = zeros([obj.grid.sz 6]);

            % aP
            coeff_p_corr(:,:,1) = -c_x(1:end-1,:) .* coeff_der_p_x(1:end-1,:,2) + c_x(2:end,:) .* coeff_der_p_x(2:end,:,1) ...
                + -c_r(:,1:end-1) .* coeff_der_p_r(:,1:end-1,2) + c_r(:,2:end) .* coeff_der_p_r(:,2:end,1);
            % aW
            coeff_p_corr(:,:,2) = c_x(1:end-1,:) .* coeff_der_p_x(1:end-1,:,1);
            % aE
            coeff_p_corr(:,:,3) = -c_x(2:end,:) .* coeff_der_p_x(2:end,:,2);
            % aS
            coeff_p_corr(:,:,4) = c_r(:,1:end-1) .* coeff_der_p_r(:,1:end-1,1);
            % aN
            coeff_p_corr(:,:,5) = -c_r(:,2:end) .* coeff_der_p_r(:,2:end,2);
            % b
            coeff_p_corr(:,:,6) = Fx(2:end,:) - Fx(1:end-1,:) + Fr(:,2:end) - Fr(:,1:end-1) ...
                + c_x(1:end-1,:) .* coeff_der_p_x(1:end-1,:,3) - c_x(2:end,:) .* coeff_der_p_x(2:end,:,3) ...
                + c_r(:,1:end-1) .* coeff_der_p_r(:,1:end-1,3) + c_r(:,2:end) .* coeff_der_p_r(:,2:end,3);
        end

        function converged = iterate(obj)
            obj.noi = obj.noi + 1;

            % Inner iterations
            for i = 1:obj.inner_iters
                % Updating properties and calculating coefficients of velocity discretized equations
                obj.update_properties();
                [coeff_vx, coeff_vr] = obj.get_coefficients_v();

                % Calculating the intermediate velocities
                vx_star=solve(coeff_vx, obj.vx, obj.solver_iters_v, randi([1 4]), obj.relaxation_factor_v);
                vr_star=solve(coeff_vr, obj.vr, obj.solver_iters_v, randi([1 4]), obj.relaxation_factor_v);

                % Calculating the face velocities using the intermediate velocities
                obj.update_face_velocities(vx_star, vr_star, coeff_vx(:,:,1), coeff_vr(:,:,1));

                % Calculating the coefficients of the pressure correction equation and solving it
                coeff_p_corr = obj.get_coefficients_p_corr(coeff_vx(:,:,1),coeff_vr(:,:,1));

                p_corr=solve(coeff_p_corr, zeros(obj.grid.sz), obj.solver_iters_p, [1;2;3;4], 1);

                % Calculating the pressure force from the pressure correction and using it to correct the velocities
                % Coefficients for pressure interpolation
                [coeff_px, coeff_pr] = linear_interpolation_scheme(obj.grid, obj.p_bds);
                % Pressure values at faces
                [px, pr] = evaluate_faces(coeff_px, coeff_pr, p_corr);
                % Multiplying the pressure by face areas
                px = px .* obj.grid.face_area_x;
                pr = pr .* obj.grid.face_area_r;
                % Adding the volume integral of pressure gradient as source terms
                p_force_x = px(1:end-1,:) - px(2:end,:);
                p_force_r = pr(:,1:end-1) - pr(:,2:end);

                vx_star = vx_star + p_force_x ./ coeff_vx(:,:,1) * obj.relaxation_factor_v;
                vr_star = vr_star + p_force_r ./ coeff_vr(:,:,1) * obj.relaxation_factor_v;

                % Updating the fields using the relaxation factor
                obj.p = obj.p + p_corr * obj.relaxation_factor_p;
                obj.vx = vx_star;
                obj.vr = vr_star;

                % Updating the face velocities and calculating the face velocities using new velocities
                obj.update_properties();
                [coeff_vx, coeff_vr] = obj.get_coefficients_v();
                obj.update_face_velocities(obj.vx, obj.vr, coeff_vx(:,:,1), coeff_vr(:,:,1));
            end

            % Updating properties
            obj.update_properties();

            % Calculating the residual and checking convergence
            [coeff_vx, coeff_vr] = obj.get_coefficients_v();
            coeff_p_corr = obj.get_coefficients_p_corr(coeff_vx(:,:,1), coeff_vr(:,:,1));

            res_vx = calculate_residual(coeff_vx, obj.vx);
            res_vr = calculate_residual(coeff_vr, obj.vr);
            res_p = norm(calculate_local_residual(coeff_p_corr, zeros(obj.grid.sz)),1) / norm(coeff_p_corr(:,:,5) .* obj.p, 1);
            
            obj.residual_history_vx(obj.noi) = res_vx;
            obj.residual_history_vr(obj.noi) = res_vr;
            obj.residual_history_p(obj.noi) = res_p;
            
            converged = (res_vx < obj.tol_v) && (res_vr < obj.tol_v) && (res_p < obj.tol_p);
        end
    end
end