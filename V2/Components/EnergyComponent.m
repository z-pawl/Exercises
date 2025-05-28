classdef EnergyComponent < IComponent
    properties
        noi (1,1)           % Number of iterations

        % Dependencies
        grid Grid2D {mustBeScalarOrEmpty}
        flow FlowComponent {mustBeScalarOrEmpty}

        temp (:,:) double   % Temperature


        % Thermophysical properties
        k (:,:) double      % Thermal conductivity
        cp (:,:) double     % Heat Capacity

        % Source terms
        q (:,:) double      % Constant heat generation term
        qt (:,:) double     % Linear heat generation term q_total = q + qt * temp

        % Functions used to calculate physical properties
        k_function (1,1) function_handle = @(x) ones(x.grid.sz);
        cp_function (1,1) function_handle = @(x) ones(x.grid.sz);

        % Functions used to calculate physical properties
        q_function (1,1) function_handle = @(x) zeros(x.grid.sz);
        qt_function (1,1) function_handle = @(x) zeros(x.grid.sz);


        % Boundary conditions
        temp_bds Boundaries {mustBeScalarOrEmpty}


        % Solver settings
        tol (1,1) double {mustBeNonnegative}            % Convergence criteria
        relaxation_factor (1,1) double                  % Relaxation factor
        inner_iters (1,1)                               % Number of inner iterations per outer iteration
        solver_iters (1,1)                              % Number of solver iterations per inner iteration
        residual_history (:,1) double                   % Residual history
    end
    methods
        function obj = EnergyComponent(grid, flow, temp, k_function, cp_function, q_function, qt_function, tol, relaxation_factor, inner_iters, solver_iters)
            arguments
                grid (1,1) Grid2D
                flow (1,1) FlowComponent
                temp (:,:) double % Initial temperature
                k_function (1,1) function_handle
                cp_function (1,1) function_handle
                q_function (1,1) function_handle
                qt_function (1,1) function_handle
                tol (1,1) double
                relaxation_factor (1,1) double
                inner_iters (1,1)
                solver_iters (1,1)
            end

            if ~isequal(size(temp), grid.sz)
                error("The initial fields must have the same size as the grid");
            end

            obj.noi = 0;
            
            obj.grid = grid;
            obj.flow = flow;

            obj.temp = temp;

            obj.k_function = k_function;
            obj.cp_function = cp_function;
            obj.q_function = q_function;
            obj.qt_function = qt_function;

            obj.temp_bds = Boundaries(grid);

            obj.tol = tol;
            obj.relaxation_factor = relaxation_factor;
            obj.inner_iters = inner_iters;
            obj.solver_iters = solver_iters;
            obj.residual_history = zeros(10000,1);
        end

        function update_properties(obj)
            obj.k = obj.k_function();
            obj.cp = obj.cp_function();
            obj.q = obj.q_function();
            obj.qt = obj.qt_function();
        end

        function coeff = get_coefficients(obj)
            % Coefficients used to calculate the physical properties at faces
            [lin_int_coeffs_x, lin_int_coeffs_r] = linear_interpolation_scheme(obj.grid, obj.grid.domain_boundary);
            
            % Properties interpolated to faces
            % Density
            [rho_x, rho_r] = evaluate_faces(lin_int_coeffs_x, lin_int_coeffs_r, obj.flow.rho);
            % Thermal conductivity
            [k_x, k_r] = evaluate_faces(lin_int_coeffs_x, lin_int_coeffs_r, obj.k);
            % Heat capacity
            [cp_x, cp_r] = evaluate_faces(lin_int_coeffs_x, lin_int_coeffs_r, obj.cp);
            

            clear lin_int_coeffs_x lin_int_coeffs_r


            % Values and normal derivatives of temperature at faces
            % Coefficients for normal derivatives
            [coeffs_x_n_der_temp, coeffs_r_n_der_temp] = central_differencing_scheme(obj.grid, obj.temp_bds);

            % Values of normal derivatives
            [n_der_temp_x, n_der_temp_r] = evaluate_faces(coeffs_x_n_der_temp, coeffs_r_n_der_temp, obj.temp);

            % Coefficients for the upwind scheme
            [coeffs_x_upwind, coeffs_r_upwind] = upwind_scheme(obj.grid, obj.flow.vx_faces, obj.flow.vr_faces, obj.temp_bds);
            
            % Coefficients for the TVD scheme
            [coeffs_x_TVD, coeffs_r_TVD] = TVD_scheme(obj.grid, obj.flow.vx_faces, obj.flow.vr_faces, @van_leer_flux_limiter, n_der_temp_x, n_der_temp_r, obj.temp_bds);

            % Coefficeints for the deferred scheme
            [coeffs_x_temp_deferred, coeffs_r_temp_deferred] = deferred_correction_face(coeffs_x_upwind, coeffs_r_upwind, coeffs_x_TVD, coeffs_r_TVD, obj.temp);


            clear n_der_temp_x n_der_temp_r coeffs_x_upwind coeffs_r_upwind coeffs_x_TVD coeffs_r_TVD;


            % F = rho * cp * v * A
            Fx = rho_x .* cp_x .* obj.flow.vx_faces .* obj.grid.face_area_x;
            Fr = rho_r .* cp_r .* obj.flow.vr_faces .* obj.grid.face_area_r;

            % D = k * A
            Dx = k_x .* obj.grid.face_area_x;
            Dr = k_r .* obj.grid.face_area_r;


            clear rho_x rho_r k_x k_r cp_x cp_r;


            coeff = assemble_coeff_array(Fx, Fr, Dx, Dr, coeffs_x_temp_deferred, coeffs_r_temp_deferred, coeffs_x_n_der_temp, coeffs_r_n_der_temp, obj.q, obj.qt, obj.grid.volume);
        end

        function converged = iterate(obj)
            obj.noi = obj.noi + 1;

            % Inner iterations
            for i = 1:obj.inner_iters
                % Updating properties and calculating coefficients of discretized equations
                obj.update_properties();
                coeff = obj.get_coefficients();
    
                % Solving the system of equations and updating the field using the relaxation factor
                obj.temp = obj.temp + obj.relaxation_factor * (solve(coeff, obj.temp, obj.solver_iters, [1; 2; 3; 4], 1) - obj.temp);
            end

            % Updating properties
            obj.update_properties();

            % Calculating the residual and checking convergence
            res = calculate_residual(obj.get_coefficients(), obj.temp);
            obj.residual_history(obj.noi) = res;
            converged = res < obj.tol;
        end
    end
end