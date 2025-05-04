classdef TemperatureCondConv2DSteady < IComponent
    properties
        sim Simulation2D {mustBeScalarOrEmpty}
        velocity IComponent {mustBeScalarOrEmpty}
        grid Grid2D {mustBeScalarOrEmpty}

        sz double               % Size

        temp double             % Temperature

        cp double               % Heat capacity
        k double                % Thermal conductivity
        q double                % Constant heat generation term
        qt double               % Linear heat generation term - q_total = q + qt*temp

        % Functions used to calculate the physical properties
        cp_function (1,1) function_handle = @k_f_default
        k_function (1,1) function_handle = @k_f_default
        q_function (1,1) function_handle = @q_f_default
        qt_function (1,1) function_handle = @qt_f_default

        boundaries_west (:,3)
        boundaries_east (:,3)
        boundaries_south (:,3)
        boundaries_north (:,3)

        tol (1,1) double {mustBePositive} = 1e-7
        residual_history (:,1) double = zeros(10000,1)
        relaxation_factor (1,1) double = 1  % Relaxation factor
        sweep_dirs (:,1) double = [1;2;3;4] % Sweep directions
        num_of_iters (1,1) double = 90      % Number of iterations per iteration
    end
    methods
            function obj = TemperatureCondConv2DSteady(sz, temp, cp_f, k_f, q_f, qt_f, tol, relaxation_factor)
                arguments
                    sz {mustBeInteger, mustBePositive}
                    temp double

                    cp_f (1,1) function_handle
                    k_f (1,1) function_handle
                    q_f (1,1) function_handle
                    qt_f (1,1) function_handle

                    tol (1,1) double {mustBePositive}
                    relaxation_factor (1,1) double {mustBePositive}
                end

                if ~isequal(sz, size(temp))
                    error("Arrays provided in the input must be the same size as the given size")
                end

                obj.sim=Simulation2D.empty;
                obj.velocity=IComponent.empty;
                obj.grid=Grid2D.empty;

                obj.sz = sz;

                obj.temp = temp;

                obj.cp_function = cp_f;
                obj.k_function = k_f;
                obj.q_function = q_f;
                obj.qt_function = qt_f;

                obj.boundaries_west=zeros(sz(2),3);
                obj.boundaries_east=zeros(sz(2),3);
                obj.boundaries_south=zeros(sz(1),3);
                obj.boundaries_north=zeros(sz(1),3);

                obj.tol = tol;
                obj.relaxation_factor = relaxation_factor;
                % obj.sweep_dirs = [1;2;3;4];
                % obj.num_of_iters = 50;

                obj.update_properties();
            end
            
            function update_properties(obj)
                obj.cp = obj.cp_function(obj);
                obj.k = obj.k_function(obj);
                obj.q = obj.q_function(obj);
                obj.qt = obj.qt_function(obj);
            end

            function coeff = get_coefficients(obj)
                % Preallocating the coefficient matrix
                coeff = zeros([obj.sz 6]);

                area_x = 2*pi*obj.grid.dr.*obj.grid.pos_cent_r;
                area_r = 2*pi*obj.grid.dx.*obj.grid.pos_stag_r;
                volume = 2*pi*obj.grid.dx.*obj.grid.dr.*obj.grid.pos_cent_r;

                % Arrays holding values of F, D and A(|P|) for every interface
                % between control volumes along x and y axis respectively
                % (for example Fx(1,1) is a value for the interface between the control volumes with indexes (1,1) and (2,1) respectively)

                % Flow rate through the cell face: F = rho*cp*v*A
                % The density and the heat capacity are lineary interpolated using the properties of neighbouring cells
                % The velocity is taken as the velocity at the cell face
                % A - area of the cell face

                Fx = Utilities.interpolate_center2faces(obj.velocity.rho,obj.grid.dx,1) ...
                    .* Utilities.interpolate_center2faces(obj.cp,obj.grid.dx,1) ...
                    .* obj.velocity.vx .* area_x;

                Fr = Utilities.interpolate_center2faces(obj.velocity.rho,obj.grid.dr,2) ...
                    .* Utilities.interpolate_center2faces(obj.cp,obj.grid.dr,2) ...
                    .* obj.velocity.vr .* area_r;


                % Conductance: D =  k/dx*A
                % dx - distance between the centers of neighbouring cells

                Dx = Utilities.interpolate_center2faces(obj.k,obj.grid.dx,1) ./ obj.grid.dx_stag .* area_x;

                Dr = Utilities.interpolate_center2faces(obj.k,obj.grid.dr,2) ./ obj.grid.dr_stag .* area_r;

                % Function A(|P|) calculated using the power law scheme
                Ax = (1-0.1*abs(Fx./Dx)).^5;
                Ar = (1-0.1*abs(Fr./Dr)).^5;

                % Every negative value is changed to 0
                Ax = Ax.*(Ax>0);
                Ar = Ar.*(Ar>0);

                % aW = Dw*A(|Pw|)+max(Fw,0)
                coeff(:,:,1) = Dx(1:end-1,:).*Ax(1:end-1,:)+Fx(1:end-1,:).*(Fx(1:end-1,:)>0);
                % aE = De*A(|Pe|)+max(-Fe,0)
                coeff(:,:,2) = Dx(2:end,:).*Ax(2:end,:)-Fx(2:end,:).*(Fx(2:end,:)<0);
                % aS = Ds*A(|Ps|)+max(Fs,0)
                coeff(:,:,3) = Dr(:,1:end-1).*Ar(:,1:end-1)+Fr(:,1:end-1).*(Fr(:,1:end-1)>0);
                % aN = Dn*A(|Pn|)+max(-Fn,0)
                coeff(:,:,4) = Dr(:,2:end).*Ar(:,2:end)-Fr(:,2:end).*(Fr(:,2:end)<0);
                % aP = aE+aW+aN+aS-SP*dx*dy
                coeff(:,:,5) = coeff(:,:,1)+coeff(:,:,2)+coeff(:,:,3)+coeff(:,:,4)-obj.qt(:,:).*volume+(Fx(2:end,:)-Fx(1:end-1,:))+(Fr(:,2:end)-Fr(:,1:end-1));
                % b = SC*dx*dy
                coeff(:,:,6) = obj.q(:,:).*volume;

                coeff = obj.apply_boundary_conditions(coeff);
            end

            function set_boundary_conditions(obj, boundary_side, boundary_index, bd_data)
                arguments
                    obj TemperatureCondConv2DSteady
                    boundary_side (1,1) string
                    boundary_index (1,1) double
                    bd_data (1,3) double
                end

                % Checking the validity of the boundary data
                if ~ismember(bd_data(1), [1 2 3])
                    disp("First value of the boundary data describes the boundary type; it has to be either 1, 2 or 3");
                    return;
                end

                % Changing the corresponding matrix holding the boundary condition
                switch boundary_side
                    case "west"
                        if 1<=boundary_index && boundary_index<=obj.sz(2)
                            obj.boundaries_west(boundary_index,:)=bd_data;
                        else
                            disp("Boundary index is out of range");
                        end
                    case "east"
                        if 1<=boundary_index && boundary_index<=obj.sz(2)
                            obj.boundaries_east(boundary_index,:)=bd_data;
                        else
                            disp("Boundary index is out of range");
                        end
                    case "south"
                        if 1<=boundary_index && boundary_index<=obj.sz(1)
                            obj.boundaries_south(boundary_index,:)=bd_data;
                        else
                            disp("Boundary index is out of range");
                        end
                    case "north"
                        if 1<=boundary_index && boundary_index<=obj.sz(1)
                            obj.boundaries_north(boundary_index,:)=bd_data;
                        else
                            disp("Boundary index is out of range");
                        end
                    otherwise
                        disp("Invalid boundary side name");
                end
            end

            function coeff = apply_boundary_conditions(obj, coeff)
                % Function calculating the changes in the coefficient for a given boundary of a cell 
                function [d_aP, d_b] = apply_boundary_condtion_cell(bd_data, aB, dq, A, k, n_dir)
                    if bd_data(1)==1
                        % Dirichlet (given temperature) boundary condition
                        d_aP=0;
                        d_b=bd_data(2)*aB;
                    elseif bd_data(1)==2
                        % Neumann (given heat flux) boundary condition
                        d_aP=-aB;
                        d_b=bd_data(2)*A*n_dir;
                    else
                        % Robin (convective boundary) boundary condition
                        m=bd_data(3)*dq+2*k;
                        d_aP=-aB*2*k/m;
                        d_b=aB*(m-2*k)/m*bd_data(2);
                    end
                end

                area_x = 2*pi*obj.grid.dr.*obj.grid.pos_cent_r;
                area_r = 2*pi*obj.grid.dx.*obj.grid.pos_stag_r;
                
                % Calculations of the boundary conditions are performed for every wall separately
                [d_aP_west, d_b_west]=arrayfun(@(x) apply_boundary_condtion_cell(obj.boundaries_west(x,:), coeff(1,x,1), obj.grid.dx(1,1), ...
                    area_x(1,x), obj.k(1,x), 1), 1:obj.sz(2));
                
                [d_aP_east, d_b_east]=arrayfun(@(x) apply_boundary_condtion_cell(obj.boundaries_east(x,:), coeff(end,x,2), obj.grid.dx(end,1), ...
                    area_x(end,x), obj.k(end,x), -1), 1:obj.sz(2));
                
                [d_aP_south, d_b_south]=arrayfun(@(x) apply_boundary_condtion_cell(obj.boundaries_south(x,:),coeff(x,1,3),obj.grid.dr(1,1), ...
                    area_r(x,1), obj.k(x,1), 1), 1:obj.sz(1));

                [d_aP_north, d_b_north]=arrayfun(@(x) apply_boundary_condtion_cell(obj.boundaries_north(x,:),coeff(x,end,4),obj.grid.dr(1,end), ...
                    area_r(x,end), obj.k(x,end), -1), 1:obj.sz(1));

                % Obtained values are used to change the coefficients of the discretized equations
                % West
                coeff(1,:,5)=coeff(1,:,5)+d_aP_west;
                coeff(1,:,6)=coeff(1,:,6)+d_b_west;
                coeff(1,:,1)=0;
                % East
                coeff(end,:,5)=coeff(end,:,5)+d_aP_east;
                coeff(end,:,6)=coeff(end,:,6)+d_b_east;
                coeff(end,:,2)=0;
                % South
                coeff(:,1,5)=coeff(:,1,5)+reshape(d_aP_south,[obj.sz(1), 1]);
                coeff(:,1,6)=coeff(:,1,6)+reshape(d_b_south, [obj.sz(1), 1]);
                coeff(:,1,3)=0;
                % North
                coeff(:,end,5)=coeff(:,end,5)+reshape(d_aP_north, [obj.sz(1), 1]);
                coeff(:,end,6)=coeff(:,end,6)+reshape(d_b_north, [obj.sz(1), 1]);
                coeff(:,end,4)=0;
            end
            
            function converged = iterate(obj)
                % Updating properties and calculating coefficients of discretized equations
                obj.update_properties();
                coeff = obj.get_coefficients();

                % Solving the system of equations and updating the temperature field using the relaxation factor
                obj.temp=obj.temp+obj.relaxation_factor*(solve(coeff,obj.temp,obj.num_of_iters,obj.sweep_dirs,1)-obj.temp);

                % Updating properties
                obj.update_properties();

                % Calculating the residual and checking convergence
                res=obj.calculate_residual(obj.get_coefficients());
                obj.residual_history(obj.sim.noi)=res;

                converged=obj.has_converged(res);
            end

            function res = calculate_residual(obj, coeff)
                % A column and a row of zeros used to pad the data in order to make sure that there is no index out of bounds error
                col0 = zeros(obj.sz(1),1);
                row0 = zeros(1,obj.sz(2));

                % r = aW*TW+aE*TE+aS*TS+aN*TN+b-aP*TP
                res_vector = [row0; coeff(2:end,:,1).*obj.temp(1:end-1,:)] + [coeff(1:end-1,:,2).*obj.temp(2:end,:); row0] ... 
                    + [col0 coeff(:,2:end,3).*obj.temp(:,1:end-1)] + [coeff(:,1:end-1,4).*obj.temp(:,2:end) col0] ...
                    + coeff(:,:,6) - coeff(:,:,5).*obj.temp;
                % s_f = || [aP1*TP1; aP2*TP2; aP3*TP3; .... aPN-1*TPN-1; aPN*TPN] ||
                scaling_factor = norm(coeff(:,:,5).*obj.temp,1);
                
                % res = || r || / s_f
                res = norm(res_vector,1)/scaling_factor;
            end

            function converged = has_converged(obj, res)
                converged = res < obj.tol;
            end

            function delete_residual_history(obj)
                obj.residual_history = zeros(10000,1);
            end
    end
end