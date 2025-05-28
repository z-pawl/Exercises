function [vx_faces, vr_faces] = rhie_chow_interpolation(grid, vx, vr, p, vx_bds, vr_bds, p_bds, aP_vx, aP_vr)
    arguments
        grid Grid2D
        vx (:,:) double
        vr (:,:) double
        p (:,:) double
        vx_bds (1,1) Boundaries
        vr_bds (1,1) Boundaries
        p_bds (1,1) Boundaries
        aP_vx (:,:) double
        aP_vr (:,:) double
    end
    %% Calculating the rhie chow correction
    % Calculating the Df
    DP_x = grid.volume ./ aP_vx;
    DP_r = grid.volume ./ aP_vr;

    % Calculating the coefficient for linear interpolation of DP and the
    % pressure gradients
    [coeffs_x, coeffs_r] = linear_interpolation_scheme(grid, grid.domain_boundary);

    % Interpolating the Df
    [Df_x, ~] = evaluate_faces(coeffs_x, coeffs_r, DP_x);
    [~, Df_r] = evaluate_faces(coeffs_x, coeffs_r, DP_r);
    clear DP_x DP_r;

    % Calculating the normal pressure derivatives at faces
    % Obtaining the coefficients
    [coeffs_x_p_der, coeffs_r_p_der] = central_differencing_scheme(grid, p_bds);
    % Calculating the derivatives
    [p_der_x, p_der_r] = evaluate_faces(coeffs_x_p_der, coeffs_r_p_der, p);
    clear coeffs_x_p_der coeffs_r_p_der;

    % Calculating the interpolated pressure derivatives at faces
    % Obtaining the coefficients for calculating pressure values at faces
    [coeffs_x_p_val, coeffs_r_p_val] = linear_interpolation_scheme(grid, p_bds);
    % Calculating the values
    [p_x, p_r] = evaluate_faces(coeffs_x_p_val, coeffs_r_p_val, p);
    clear coeffs_x_p_val coeffs_r_p_val;
    % Calculating the derivatives at cell centers
    p_der_center_x = grid.face_area_x(2:end,:) .* p_x(2:end,:) - grid.face_area_x(1:end-1,:) .* p_x(1:end-1,:);
    p_der_center_r = grid.face_area_r(:,2:end) .* p_r(:,2:end) - grid.face_area_r(:,1:end-1) .* p_r(:,1:end-1);
    clear p_x p_r;
    % Calculating the interpolated derivatives
    [p_der_x_int, ~] = evaluate_faces(coeffs_x, coeffs_r, p_der_center_x);
    [~, p_der_r_int] = evaluate_faces(coeffs_x, coeffs_r, p_der_center_r);
    clear p_der_center_x p_der_center_r;

    
    % Calculating the correction
    correction_x = Df_x .* (p_der_x_int - p_der_x);
    correction_r = Df_r .* (p_der_r_int - p_der_r);
    clear p_der_x p_der_r p_der_x_int p_der_r_int;

    %% Obtaining the interpolation
    % Calculating the coefficient for linear interpolation of velocities
    [coeffs_x, ~] = linear_interpolation_scheme(grid, vx_bds);
    [~, coeffs_r] = linear_interpolation_scheme(grid, vr_bds);

    % Applying the correction
    coeffs_x(:,:,3) = coeffs_x(:,:,3) + correction_x;
    coeffs_r(:,:,3) = coeffs_r(:,:,3) + correction_r;
    clear correction_x correction_r;
    
    % Applying the boundary conditions
    [coeffs_x, ~] = vx_bds.apply_boundary_condition_value(coeffs_x, coeffs_r);
    [~, coeffs_r] = vr_bds.apply_boundary_condition_value(coeffs_x, coeffs_r);

    %% Calculating face values
    [vx_faces, ~] = evaluate_faces(coeffs_x, coeffs_r, vx);
    [~, vr_faces] = evaluate_faces(coeffs_x, coeffs_r, vr);
end