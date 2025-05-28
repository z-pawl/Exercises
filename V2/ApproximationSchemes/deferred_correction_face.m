function [coeffs_x, coeffs_r] = deferred_correction_face(coeffs_x_low, coeffs_r_low, coeffs_x_high, coeffs_r_high, field)
    arguments
        coeffs_x_low (:,:,3) double
        coeffs_r_low (:,:,3) double
        coeffs_x_high (:,:,3) double
        coeffs_r_high (:,:,3) double
        field (:,:) double
    end
    coeffs_x = coeffs_x_low;
    coeffs_r = coeffs_r_low;

    [delta_vals_x, delta_vals_r] = evaluate_faces(coeffs_x_high - coeffs_x_low, coeffs_r_high - coeffs_r_low, field);


    coeffs_x(:,:,3) = coeffs_x(:,:,3) + delta_vals_x;
    coeffs_r(:,:,3) = coeffs_r(:,:,3) + delta_vals_r;
end