classdef Grid2D < handle
    properties
        % Size of the grid
        sz (1,2) double

        % Cell sizes along given direction
        dx (:,1) double
        dr (1,:) double

        % Cell center positions along given direction
        cent_pos_x (:,1) double
        cent_pos_r (1,:) double

        % Position of the center of faces perpendicular to given direction
        face_pos_x (:,1) double
        face_pos_r (1,:) double

        % Area of faces perpendicular to given direction
        face_area_x (:,:) double
        face_area_r (:,:) double

        % Cell volumes
        volume (:,:) double

        % Domain boundary
        domain_boundary Boundaries {mustBeScalarOrEmpty}
    end
    methods
        function obj = Grid2D(dx, dr, offset_r)
            arguments
                dx (:,1) double % Cell sizes along the x direction
                dr (1,:) double % Cell sizes along the r direction
                offset_r (1,1) double % Offset of the first cell face along the r direction
            end
            % Assigning the properties
            obj.sz = [size(dx,1) size(dr,2)];

            obj.dx = dx;
            obj.dr = dr;

            obj.face_pos_x = cumsum([0; dx], 1);
            obj.face_pos_r = cumsum([offset_r dr], 2);

            obj.cent_pos_x = (obj.face_pos_x(1:end-1,:) + obj.face_pos_x(2:end,:)) / 2;
            obj.cent_pos_r = (obj.face_pos_r(:,1:end-1) + obj.face_pos_r(:,2:end)) / 2;

            obj.face_area_x = repmat(2 * pi .* dr .* obj.cent_pos_r, obj.sz(1) + 1, 1);
            obj.face_area_r = 2 * pi .* dx .* obj.face_pos_r;

            obj.volume = 2 * pi .* dx .* dr .* obj.cent_pos_r;

            obj.domain_boundary = Boundaries(obj);
        end
    end
end