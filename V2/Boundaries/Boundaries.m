classdef Boundaries < handle
    properties
        grid Grid2D {mustBeScalarOrEmpty}

        % Boundary conditions for faces normal to a specific coordinate direction.
        % Each face is described by 5 variables:
        % is_bd, orientation, a, b, c
        % - is_bd: whether given face is a boundary
        % - orientation: 1 if the inward normal vector aligns with the coordinate direction, -1 otherwise
        % - a, b, c: values specifying the Robin bd condition a*val_f + b*val'_f = c (other bd conditions are special cases)

        boundary_conditions_x struct
        boundary_conditions_r struct
    end
    methods
        function obj = Boundaries(grid)
            arguments
                grid Grid2D
            end

            obj.grid = grid;

            sz_x = grid.sz + [1 0];
            sz_r = grid.sz + [0 1];

            obj.boundary_conditions_x = struct("is_bd", false(sz_x), ...
                "orientation", zeros(sz_x), ...
                "a", zeros(sz_x), ...
                "b", zeros(sz_x), ...
                "c", zeros(sz_x));
            
            obj.boundary_conditions_r = struct("is_bd", false(sz_r), ...
                "orientation", zeros(sz_r), ...
                "a", zeros(sz_r), ...
                "b", zeros(sz_r), ...
                "c", zeros(sz_r));
        end
        function set_boundary_condition(obj, direction, face_index, orientation, a, b, c)
            arguments
                obj Boundaries
                direction (1,1) string
                face_index (1,2) double
                orientation (1,1) double
                a (1,1) double
                b (1,1) double
                c (1,1) double
            end

            if orientation ~= -1 && orientation ~= 1
                error("The orientation has to be either -1 or 1");
            end

            switch direction
                case "x"
                    if (1 > face_index(1) || face_index(1) > obj.grid.sz(1)+1) || (1 > face_index(2) || face_index(2) > obj.grid.sz(2))
                        error("Provided face index is out of bounds");
                    end
                    obj.boundary_conditions_x.is_bd(face_index(1),face_index(2)) = true;
                    obj.boundary_conditions_x.orientation(face_index(1),face_index(2)) = orientation;
                    obj.boundary_conditions_x.a(face_index(1),face_index(2)) = a;
                    obj.boundary_conditions_x.b(face_index(1),face_index(2)) = b;
                    obj.boundary_conditions_x.c(face_index(1),face_index(2)) = c;
                case "r"
                    if (1 > face_index(1) || face_index(1) > obj.grid.sz(1)) || (1 > face_index(2) || face_index(2) > obj.grid.sz(2)+1)
                        error("Provided face index is out of bounds");
                    end
                    obj.boundary_conditions_r.is_bd(face_index(1),face_index(2)) = true;
                    obj.boundary_conditions_r.orientation(face_index(1),face_index(2)) = orientation;
                    obj.boundary_conditions_r.a(face_index(1),face_index(2)) = a;
                    obj.boundary_conditions_r.b(face_index(1),face_index(2)) = b;
                    obj.boundary_conditions_r.c(face_index(1),face_index(2)) = c;
                otherwise
                    error("Provided direction has to be a name of one of the directions, either 'x' or 'r'");
            end
        end
        
        % Function applying the boundary conditions to field values at boundaries
        function [coeffs_x, coeffs_r] = apply_boundary_condition_value(obj, coeffs_x, coeffs_r)
            arguments
                obj Boundaries
                % Coefficient representing the boundary value in the form
                % [aP, aN, b] where P means the previous cell and N the
                % next one The value at boundary is constructed as follows:
                % val_face = aP*phiP + aN*phiN + b
                coeffs_x (:,:,3) double
                coeffs_r (:,:,3) double
            end

            % Applying the bd conditions to the x faces
            for i = 1:obj.grid.sz(1)+1
                for j = 1:obj.grid.sz(2)
                    if obj.boundary_conditions_x.is_bd(i,j)
                        orientation = obj.boundary_conditions_x.orientation(i,j);
                        a = obj.boundary_conditions_x.a(i,j);
                        b = obj.boundary_conditions_x.b(i,j);
                        c = obj.boundary_conditions_x.c(i,j);

                        if orientation == 1
                            dx = obj.grid.dx(i,1) / 2;
                            denominator = a * dx * orientation - b;

                            coeffs_x(i,j,1) = 0;
                            coeffs_x(i,j,2) = -b / denominator;
                            coeffs_x(i,j,3) = c * dx * orientation / denominator;
                        else
                            dx = obj.grid.dx(i-1,1) / 2;
                            denominator = a * dx * orientation - b;

                            coeffs_x(i,j,1) = -b / denominator;
                            coeffs_x(i,j,2) = 0;
                            coeffs_x(i,j,3) = c * dx * orientation / denominator;
                        end
                    end
                end
            end

            % Applying the bd conditions to the r faces
            for i = 1:obj.grid.sz(1)
                for j = 1:obj.grid.sz(2)+1
                    if obj.boundary_conditions_r.is_bd(i,j)
                        orientation = obj.boundary_conditions_r.orientation(i,j);
                        a = obj.boundary_conditions_r.a(i,j);
                        b = obj.boundary_conditions_r.b(i,j);
                        c = obj.boundary_conditions_r.c(i,j);

                        if orientation == 1
                            dr = obj.grid.dr(1,j) / 2;
                            denominator = a * dr * orientation - b;

                            coeffs_r(i,j,1) = 0;
                            coeffs_r(i,j,2) = -b / denominator;
                            coeffs_r(i,j,3) = c * dr * orientation / denominator;
                        else
                            dr = obj.grid.dr(1,j-1) / 2;
                            denominator = a * dr * orientation - b;

                            coeffs_r(i,j,1) = -b / denominator;
                            coeffs_r(i,j,2) = 0;
                            coeffs_r(i,j,3) = c * dr * orientation / denominator;
                        end
                    end
                end
            end
        end

        % Function applying the boundary conditions to field normal derivatives at boundaries
        function [coeffs_x, coeffs_r] = apply_boundary_condition_normal_derivative(obj, coeffs_x, coeffs_r)
            arguments
                obj Boundaries
                % Coefficient representing the boundary derivative in the form
                % [aP, aN, b] where P means the previous cell and N the
                % next one The value at boundary is constructed as follows:
                % derivative_face = aP*phiP + aN*phiN + b
                coeffs_x (:,:,3) double
                coeffs_r (:,:,3) double
            end

            % Applying the bd conditions to the x faces
            for i = 1:obj.grid.sz(1)+1
                for j = 1:obj.grid.sz(2)
                    if obj.boundary_conditions_x.is_bd(i,j)
                        orientation = obj.boundary_conditions_x.orientation(i,j);
                        a = obj.boundary_conditions_x.a(i,j);
                        b = obj.boundary_conditions_x.b(i,j);
                        c = obj.boundary_conditions_x.c(i,j);

                        if orientation == 1
                            dx = obj.grid.dx(i,1) / 2;
                            denominator = a * dx * orientation - b;

                            coeffs_x(i,j,1) = 0;
                            coeffs_x(i,j,2) = a / denominator;
                            coeffs_x(i,j,3) = c / denominator;
                        else
                            dx = obj.grid.dx(i-1,1) / 2;
                            denominator = a * dx * orientation - b;

                            coeffs_x(i,j,1) = a / denominator;
                            coeffs_x(i,j,2) = 0;
                            coeffs_x(i,j,3) = c / denominator;
                        end
                    end
                end
            end

            % Applying the bd conditions to the r faces
            for i = 1:obj.grid.sz(1)
                for j = 1:obj.grid.sz(2)+1
                    if obj.boundary_conditions_r.is_bd(i,j)
                        orientation = obj.boundary_conditions_r.orientation(i,j);
                        a = obj.boundary_conditions_r.a(i,j);
                        b = obj.boundary_conditions_r.b(i,j);
                        c = obj.boundary_conditions_r.c(i,j);

                        if orientation == 1
                            dr = obj.grid.dr(1,j) / 2;
                            denominator = a * dr * orientation - b;

                            coeffs_r(i,j,1) = 0;
                            coeffs_r(i,j,2) = a / denominator;
                            coeffs_r(i,j,3) = c / denominator;
                        else
                            dr = obj.grid.dr(1,j-1) / 2;
                            denominator = a * dr * orientation - b;

                            coeffs_r(i,j,1) = a / denominator;
                            coeffs_r(i,j,2) = 0;
                            coeffs_r(i,j,3) = c / denominator;
                        end
                    end
                end
            end
        end
    end
end