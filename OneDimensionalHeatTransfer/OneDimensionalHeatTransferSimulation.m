% Simulates one-dimensional conduction and convection for given velocity field using
% discretization equations described in "Patankar, S. Numerical Heat Transfer and Fluid Flow"

classdef OneDimensionalHeatTransferSimulation
    properties
        % Number of control volumes in the simulation
        sim_size (1,1) {mustBeInteger, mustBeNonnegative}
        
        % Vectors of properties for each of the control volumes
        dx (:,1) double             % Length of the control volume
        k (:,1) double              % Thermal conductivity
        dcp (:,1) double            % Product of density and specific heat capacity
        q (:,1) double              % Constant heat generation term
        qt (:,1) double             % Linear heat generation term - q_total = q + qt*temp
        temp (:,1) double           % Temperature

        % A collocated grid approach is used to represent the velocity field
        % throughout the simulation domain
        % No tests are performed to check if the velocity field satisfies
        % the continuity equation
        u (:,1) double              % Velocity

        % Functions used to calculate the physical properties
        k_function (1,1) function_handle = @k_f_default
        dcp_function (1,1) function_handle = @dcp_f_default
        q_function (1,1) function_handle = @q_f_default
        qt_function (1,1) function_handle = @qt_f_default

        % Coefficients for the boundary conditions discretization equations in the form:
        % aB*TB = aI*TI + b
        % aB, TB - coefficient and temperature of the boundary cell
        % aI, TI - coefficient and temperature of the neighbouring cell
        % b - constant
        % Boundaries are represented as [aB; aI; b]
        left_boundary (3,1) double = [0; 1; 0]
        right_boundary (3,1) double = [0; 1; 0]

        noi (1,1) {mustBeInteger, mustBeNonnegative}    % Number of performed iterations

        tol (1,1) double {mustBePositive} = 1e-7                % Error tolerance
        max_iters (1,1) {mustBeInteger, mustBePositive} = 100   % The aximum number of iterations for given calculation
        rel_change_history (:,1) double = zeros(10000,1)               % Vector storing the values of the errors for previous iterations
    end
    methods
        % Constructor method
        function obj = OneDimensionalHeatConductionSimulation(sim_size, dx, temp, dt, tol, max_iters, ...
                k_function, dcp_function, q_function, qt_function)
            arguments (Input)
                sim_size (1,1) {mustBeInteger, mustBeNonnegative}
                dx (:,1) double
                temp (:,1) double

                dt (1,1) double

                tol (1,1) double {mustBePositive}
                max_iters (1,1) {mustBeInteger, mustBePositive}

                k_function (1,1) function_handle
                dcp_function (1,1) function_handle
                q_function (1,1) function_handle
                qt_function (1,1) function_handle
            end
            arguments (Output)
                obj OneDimensionalHeatConductionSimulation
            end
            
            obj.sim_size = sim_size;
            obj.dx = dx;
            obj.temp = temp;

            obj.noi = 0;
            obj.dt = dt;
            obj.elapsed = 0;

            obj.tol = tol;
            obj.max_iters = max_iters;

            obj.k_function = k_function;
            obj.dcp_function = dcp_function;
            obj.q_function = q_function;
            obj.qt_function = qt_function;

            obj = obj.update_properties();
        end

        % Calculates the position of each node
        function pos = get_nodes_position(obj)
            pos = zeros(obj.sim_size,1);
            pos(1) = obj.dx(1)/2; 
            for i = 2:obj.sim_size
                pos(i) = pos(i-1)+(obj.dx(i-1)+obj.dx(i))/2;
            end
        end
        
        % Calculates the coefficients of the discretization equations for
        % each node for the steady state
        % Discretization equations are of the form:
        % aP*TP = aW*TW + aE*TE + b
        % aP, TP - coefficient and temperature of given control volume
        % aW, TW - coefficient and temperature of the control volume west to the given one
        % aE, TE - coefficient and temperature of the control volume east to the given one
        % b - constant
        % Coefficient matrix holds those values as colums in the following
        % order: [aW aE aP b]
        function coeff = get_coefficients_steady(obj)
            % Preallocating the matrix
            coeff = zeros(obj.sim_size, 4);
            % aW
            coeff(2:end-1,1) = 2*(obj.dx(1:end-2)./obj.k(1:end-2)+obj.dx(2:end-1)./obj.k(2:end-1)).^(-1);
            % aE
            coeff(2:end-1,2) = 2*(obj.dx(2:end-1)./obj.k(2:end-1)+obj.dx(3:end)./obj.k(3:end)).^(-1);
            % aP
            coeff(2:end-1,3) = coeff(2:end-1,1)+coeff(2:end-1,2)-obj.qt(2:end-1).*obj.dx(2:end-1);
            % b
            coeff(2:end-1,4) = obj.q(2:end-1).*obj.dx(2:end-1);

            % Left boundary
            coeff(1,2) = obj.left_boundary(1);
            coeff(1,3) = obj.left_boundary(2);
            coeff(1,4) = obj.left_boundary(3);

            % Right boundary
            coeff(end,1) = obj.right_boundary(1);
            coeff(end,3) = obj.right_boundary(2);
            coeff(end,4) = obj.right_boundary(3);
        end

        % Calculates the coefficients of the discretization equations for
        % each node for the transient state
        function coeff = get_coefficients_unsteady(obj, old_temps, old_dcp)
            % Preallocating the matrix
            coeff = zeros(obj.sim_size, 4);
            % aW
            coeff(2:end-1,1) = 2*(obj.dx(1:end-2)./obj.k(1:end-2)+obj.dx(2:end-1)./obj.k(2:end-1)).^(-1);
            % aE
            coeff(2:end-1,2) = 2*(obj.dx(2:end-1)./obj.k(2:end-1)+obj.dx(3:end)./obj.k(3:end)).^(-1);
            % b
            coeff(2:end-1,4) = obj.q(2:end-1).*obj.dx(2:end-1)+old_dcp(2:end-1).*obj.dx(2:end-1)/obj.dt.*old_temps(2:end-1);
            % aP
            coeff(2:end-1,3) = coeff(2:end-1,1)+coeff(2:end-1,2)+obj.dcp(2:end-1).*obj.dx(2:end-1)/obj.dt-obj.qt(2:end-1).*obj.dx(2:end-1);

            % Left boundary
            coeff(1,2) = obj.left_boundary(1);
            coeff(1,3) = obj.left_boundary(2);
            coeff(1,4) = obj.left_boundary(3);

            % Right boundary
            coeff(end,1) = obj.right_boundary(1);
            coeff(end,3) = obj.right_boundary(2);
            coeff(end,4) = obj.right_boundary(3);
        end
        
        % Calculates the relative change between the old and new
        % temperature. 
        function rel_change = calculate_relative_change(obj, new_temp)
            rel_change = norm(obj.temp-new_temp)/norm(obj.temp);
        end

        % Updates each of the following properties:
        function obj = update_properties(obj)
            obj = obj.k_function(obj);
            obj = obj.dcp_function(obj);
            obj = obj.q_function(obj);
            obj = obj.qt_function(obj);
        end

        % Performs the calculations for the steady state
        function obj = steady(obj)
            % Number of iterations before calculations
            noi_start=obj.noi;

            % Calculations are performed while the number of iterations is
            % smaller than the limit
            while obj.noi-noi_start<obj.max_iters

                % Incrementing the number of iterations
                obj.noi=obj.noi+1;

                % Recalculating the physical properties
                obj = obj.update_properties();

                % Assigning the coefficients to the vectors for the solver
                coeff = obj.get_coefficients_steady();
                a = -coeff(2:end,1);
                b = coeff(1:end,3);
                c = -coeff(1:end-1,2);
                d = coeff(1:end,4);
                
                % Calculating the solution
                new_temp = tdma(a,b,c,d);
    
                % Calculating relative change and checking convergence
                obj.rel_change_history(obj.noi) = obj.calculate_relative_change(new_temp);
                converged = obj.rel_change_history(obj.noi)<obj.tol;

                % Updating the temperature
                obj.temp = new_temp;

                if converged
                    break
                end
            end
        end
        
        % Performs the calculations for the transient state
        function obj = unsteady(obj, timesteps)
            arguments
                obj 
                timesteps (1,1) {mustBeInteger,mustBePositive} 
            end

            % Calculations are performed for the given amount of timesteps
            for i = 1:timesteps

                % Incrementing the elapsed time
                
                obj.elapsed = obj.elapsed+obj.dt;

                % Number of iterations before calculations
                noi_start=obj.noi;

                % Old values of dcp and temperature
                old_dcp = obj.dcp;
                old_temp = obj.temp;

                % Calculations are performed while the number of iterations is
                % smaller than the limit
                while obj.noi-noi_start<obj.max_iters

                    % Incrementing the number of iterations
                    obj.noi = obj.noi+1;

                    % Recalculating the physical properties
                    obj = obj.update_properties();

                    % Assigning the coefficients to the vectors for the solver
                    coeff = obj.get_coefficients_unsteady(old_temp, old_dcp);
                    a = -coeff(2:end,1);
                    b = coeff(1:end,3);
                    c = -coeff(1:end-1,2);
                    d = coeff(1:end,4);
                    
                    % Calculating the solution
                    new_temp = tdma(a,b,c,d);
        
                    % Calculating relative change and checking convergence
                    obj.rel_change_history(obj.noi) = obj.calculate_relative_change(new_temp);
                    converged = obj.rel_change_history(obj.noi)<obj.tol;
    
                    % Updating the temperature
                    obj.temp = new_temp;
    
                    if converged
                        break
                    end
                end
            end
        end

        % Displays the relative change history
        function plot_rel_change_history(obj)
            semilogy(1:obj.noi,obj.rel_change_history(1:obj.noi)+1e-20);
            ylabel("Relative change");
            xlabel("Iteration");
        end

        % Deletes the relative change history
        function obj = delete_rel_change_history(obj)
            obj.noi = 0;
            obj.elapsed = 0;
            obj.rel_change_history = zeros(10000,1);
        end
    
        % Displays the temperature
        function plot_temperature(obj)
            % Calculationg the positions of the nodes
            x = get_nodes_position(obj);
            
            % Plotting the graph
            plot(x, obj.temp);
            xlabel("Position [m]");
            ylabel("Temperature [K]");
        end
    end

    % Functions changing the boundary conditions
    % boundary_side describes which boundary to change
    % Value of 0 refers to the left boundary and the value of 1 refers to
    % the right one
    methods
        function obj = set_dirichlet_boundary_condition(obj, boundary_side, temperature)
            arguments
                obj 
                boundary_side (1,1) {mustBeInteger}
                temperature (1,1) double
            end
            if boundary_side == 0
                obj.left_boundary = [0; 1; temperature];
            elseif boundary_side == 1
                obj.right_boundary = [0; 1; temperature];
            else
                error(['The boundary argument has to be either 0 or 1', newline, ...
                    '0 represents the left boundary and 1 represents the right one']);
            end
        end

        function obj = set_neumann_boundary_condition(obj, boundary_side, heat_flux)
            arguments
                obj 
                boundary_side (1,1) {mustBeInteger}
                heat_flux (1,1) double
            end
            if boundary_side == 0
                aI = 2/(obj.dx(1)/obj.k(1)+obj.dx(2)/obj.k(2));
                aB = aI-obj.qt(1)*obj.dx(1);
                b = obj.q(1)*obj.dx(1)+heat_flux;
                obj.left_boundary = [aI; aB; b];
            elseif boundary_side == 1
                aI = 2/(obj.dx(1)/obj.k(1)+obj.dx(2)/obj.k(2));
                aB = aI-obj.qt(1)*obj.dx(1);
                b = obj.q(1)*obj.dx(1)-heat_flux;
                obj.right_boundary = [aI; aB; b];
            else
                error(['The boundary argument has to be either 0 or 1', newline, ...
                    '0 represents the left boundary and 1 represents the right one']);
            end
        end

        function obj = set_robin_boundary_condition(obj, boundary_side, heat_transfer_coefficient, fluid_temperature)
            arguments
                obj 
                boundary_side (1,1) {mustBeInteger}
                heat_transfer_coefficient (1,1) double
                fluid_temperature (1,1) double
            end
            if boundary_side == 0
                aI = 2/(obj.dx(1)/obj.k(1)+obj.dx(2)/obj.k(2));
                aB = aI-obj.qt(1)*obj.dx(1)+heat_transfer_coefficient;
                b = obj.q(1)*obj.dx(1)+heat_transfer_coefficient*fluid_temperature;
                obj.left_boundary = [aI; aB; b];
            elseif boundary_side == 1
                aI = 2/(obj.dx(1)/obj.k(1)+obj.dx(2)/obj.k(2));
                aB = aI-obj.qt(1)*obj.dx(1)-heat_transfer_coefficient;
                b = obj.q(1)*obj.dx(1)-heat_transfer_coefficient*fluid_temperature;
                obj.right_boundary = [aI; aB; b];
            else
                error(['The boundary argument has to be either 0 or 1', newline, ...
                    '0 represents the left boundary and 1 represents the right one']);
            end
        end
    end
end

% Default functions
% k_f_default and dcp_f_default - value equal to one throughout the whole simulation domain
% q_f_default and q_t_default - value equal to zero throughout the whole simulation domain
function obj = k_f_default(obj)
    arguments
        obj OneDimensionalHeatConductionSimulation
    end
    obj.k = ones(obj.sim_size,1);
end

function obj = dcp_f_default(obj)
    arguments
        obj OneDimensionalHeatConductionSimulation
    end
    obj.dcp = ones(obj.sim_size,1);
end

function obj = q_f_default(obj)
    arguments
        obj OneDimensionalHeatConductionSimulation
    end
    obj.q = zeros(obj.sim_size,1);
end

function obj = qt_f_default(obj)
    arguments
        obj OneDimensionalHeatConductionSimulation
    end
    obj.qt = zeros(obj.sim_size,1);
end