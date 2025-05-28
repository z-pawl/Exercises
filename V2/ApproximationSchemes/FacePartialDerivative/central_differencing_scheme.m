% Function returning coefficients used to compute normal derivatives at
% cell faces.
function [coeffs_x, coeffs_r] = central_differencing_scheme(grid, field_bd_conditions)
    arguments
        grid (1,1) Grid2D
        field_bd_conditions (1,1) Boundaries
    end
    sz = grid.sz;

    coeffs_x = zeros([(sz + [1 0]) 3]);
    coeffs_r = zeros([(sz + [0 1]) 3]);

    % Inverse of the distance between cell centers along the x direction
    center_inv_distance_x = repmat(2 ./ ([0; grid.dx] + [grid.dx; 0]), 1, sz(2));

    coeffs_x(:,:,1) = -center_inv_distance_x;
    coeffs_x(:,:,2) = center_inv_distance_x;

    % Inverse of the distance between cell centers along the r direction
    center_inv_distance_r = repmat(2 ./ ([0 grid.dr] + [grid.dr 0]), sz(1), 1);

    coeffs_r(:,:,1) = -center_inv_distance_r;
    coeffs_r(:,:,2) = center_inv_distance_r;

    % Applying the boundary conditions
    [coeffs_x, coeffs_r] = field_bd_conditions.apply_boundary_condition_normal_derivative(coeffs_x, coeffs_r);
end