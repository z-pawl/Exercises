classdef Grid2D < handle
    properties
        % Size of the grid
        sz 

        % Size of a cell along given dimension
        dx (:,1) double
        dr (1,:) double

        % Size of a staggered cell staggered along given dimension
        dx_stag (:,1) double
        dr_stag (1,:) double

        % Position of the center of the cell w.r.t. the outer vertex of the (1,1) cell
        pos_cent_x (:,1) double
        pos_cent_r (1,:) double

        % Position of the center of the grid staggered along the x and y direction respectively w.r.t. the outer vertex of the (1,1) cell
        pos_stag_x (:,1) double
        pos_stag_r (1,:) double
    end
    methods
        function obj = Grid2D(dx, dr, offset_r)
            arguments
                dx (:,1) double
                dr (1,:) double
                offset_r (1,1) double
            end
            % Assigning the properties
            obj.sz = [size(dx,1) size(dr,2)];

            obj.dx = dx;
            obj.dr = dr;

            obj.dx_stag = ([0; dx]+[dx; 0])/2;
            obj.dr_stag = ([0 dr]+[dr 0])/2;

            obj.pos_stag_x = cumsum([0;dx]);
            obj.pos_stag_r = cumsum([offset_r dr]);
            
            obj.pos_cent_x = (obj.pos_stag_x(1:end-1,1)+obj.pos_stag_x(2:end,1))/2;
            obj.pos_cent_r = (obj.pos_stag_r(1,1:end-1)+obj.pos_stag_r(1,2:end))/2;
        end
    end
end