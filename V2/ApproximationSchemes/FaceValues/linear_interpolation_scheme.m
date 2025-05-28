function [coeffs_x, coeffs_r] = linear_interpolation_scheme(grid, field_bd_conditions)
    arguments
        grid (1,1) Grid2D
        field_bd_conditions (1,1) Boundaries 
    end
    sz = grid.sz;

    coeffs_x = zeros([(sz + [1 0]) 3]);
    coeffs_r = zeros([(sz + [0 1]) 3]);

    % Inverse of twice the distance between cell centers along the x direction
    center_inv_distance_x = repmat(1 ./ ([1e30; grid.dx] + [grid.dx; 1e30]), 1, sz(2));

    coeffs_x(:,:,1) = [grid.dx; 1e30] .* center_inv_distance_x;
    coeffs_x(:,:,2) = [1e30; grid.dx] .* center_inv_distance_x;

    % Inverse of twice the distance between cell centers along the r direction
    center_inv_distance_r = repmat(1 ./ ([1e30 grid.dr] + [grid.dr 1e30]), sz(1), 1);

    coeffs_r(:,:,1) = [grid.dr 1e30] .* center_inv_distance_r;
    coeffs_r(:,:,2) = [1e30 grid.dr] .* center_inv_distance_r;

    % Applying the boundary conditions
    [coeffs_x, coeffs_r] = field_bd_conditions.apply_boundary_condition_value(coeffs_x, coeffs_r);
end