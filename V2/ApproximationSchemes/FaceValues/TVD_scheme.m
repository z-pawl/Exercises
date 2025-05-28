function [coeffs_x, coeffs_r] = TVD_scheme(grid, vel_x_faces, vel_r_faces, flux_limiter, field_face_derivatives_x, field_face_derivatives_r, field_bd_conditions)
    arguments
        grid (1,1) Grid2D
        vel_x_faces (:,:) double
        vel_r_faces (:,:) double
        flux_limiter (1,1) function_handle
        field_face_derivatives_x (:,:) double 
        field_face_derivatives_r (:,:) double
        field_bd_conditions (1,1) Boundaries
    end
    sz = grid.sz;

    coeffs_x = zeros([(sz + [1 0]) 3]);
    coeffs_r = zeros([(sz + [0 1]) 3]);

    for i = 2:sz(1)
        for j = 1:sz(2)
            if vel_x_faces(i,j) >= 0
                limiter = flux_limiter(field_face_derivatives_x(i,j), field_face_derivatives_x(i-1,j));
                coeffs_x(i,j,1) = 1 - 0.5 * limiter;
                coeffs_x(i,j,2) = 0.5 * limiter;
                % coeffs_x(i,j,3) = 0;
            else
                limiter = flux_limiter(field_face_derivatives_x(i,j), field_face_derivatives_x(i+1,j));
                coeffs_x(i,j,1) = 0.5 * limiter;
                coeffs_x(i,j,2) = 1 - 0.5 * limiter;
                % coeffs_x(i,j,3) = 0;
            end
        end
    end

    for i = 1:sz(1)
        for j = 2:sz(2)
            if vel_r_faces(i,j) >= 0
                limiter = flux_limiter(field_face_derivatives_r(i,j), field_face_derivatives_r(i,j-1));
                coeffs_r(i,j,1) = 1 - 0.5 * limiter;
                coeffs_r(i,j,2) = 0.5 * limiter;
                % coeffs_x(i,j,3) = 0;
            else
                limiter = flux_limiter(field_face_derivatives_r(i,j), field_face_derivatives_r(i,j+1));
                coeffs_r(i,j,1) = 0.5 * limiter;
                coeffs_r(i,j,2) = 1 - 0.5 * limiter;
                % coeffs_x(i,j,3) = 0;
            end
        end
    end

    [coeffs_x, coeffs_r] = field_bd_conditions.apply_boundary_condition_value(coeffs_x, coeffs_r);
end