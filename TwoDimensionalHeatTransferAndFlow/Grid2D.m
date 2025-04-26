classdef Grid2D < handle
    properties
        % Size of the grid
        sz 

        % Size of a cell along given dimension
        dx (:,1) double
        dy (1,:) double

        % Size of a staggered cell staggered along given dimension
        dx_stag (:,1) double
        dy_stag (1,:) double

        % Position of the center of the cell w.r.t. the outer vertex of the (1,1) cell
        pos_cent_x (:,1) double
        pos_cent_y (1,:) double

        % Position of the center of the grid staggered along the x and y direction respectively w.r.t. the outer vertex of the (1,1) cell
        pos_stag_x (:,1) double
        pos_stag_y (1,:) double
    end
    methods
        function obj = Grid2D(dx, dy)
            arguments
                dx (:,1) double
                dy (1,:) double
            end
            % Assigning the properties
            obj.sz = [size(dx,1) size(dy,2)];

            obj.dx = dx;
            obj.dy = dy;

            obj.dx_stag = ([0; dx]+[dx; 0])/2;
            obj.dy_stag = ([0 dy]+[dy 0])/2;

            obj.pos_stag_x = cumsum([0;dx]);
            obj.pos_stag_y = cumsum([0 dy]);
            
            obj.pos_cent_x = (obj.pos_stag_x(1:end-1,1)+obj.pos_stag_x(2:end,1))/2;
            obj.pos_cent_y = (obj.pos_stag_y(1,1:end-1)+obj.pos_stag_y(1,2:end))/2;
        end
    end
end