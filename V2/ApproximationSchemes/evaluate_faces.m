% Function calculates values at faces using the provided coefficients [aP, aN, b]
% val_f = aP * val_P + aN * val_N + b
function [vals_x, vals_r] = evaluate_faces(coeffs_x, coeffs_r, field)
    arguments
        coeffs_x (:,:,3) double
        coeffs_r (:,:,3) double
        field (:,:) double
    end
    sz = size(field);

    % A column and a row of zeros used to pad the data in order to make sure that there is no index out of bounds error
    col0 = zeros(sz(1),1);
    row0 = zeros(1,sz(2));

    vals_x = coeffs_x(:,:,1) .* [row0; field] + coeffs_x(:,:,2) .* [field; row0] + coeffs_x(:,:,3);
    vals_r = coeffs_r(:,:,1) .* [col0 field] + coeffs_r(:,:,2) .* [field col0] + coeffs_r(:,:,3);
end