function coeff = assemble_coeff_array(Fx, Fr, Dx, Dr, coeffs_value_x, coeffs_value_r, coeffs_derivative_x, coeffs_derivative_r, src, src_lin, volume)
    arguments
        Fx (:,:) double
        Fr (:,:) double
        Dx (:,:) double
        Dr (:,:) double
        coeffs_value_x (:,:,3) double
        coeffs_value_r (:,:,3) double
        coeffs_derivative_x (:,:,3) double
        coeffs_derivative_r (:,:,3) double
        src (:,:) double
        src_lin (:,:) double
        volume (:,:) double
    end

    % Preallocating the array
    % Coefficients are stored in the following order:
    % [aP aW aE aS aN b]
    coeff = zeros([size(volume) 6]);

    % aP
    coeff(:,:,1)  = Fx(2:end,:) .* coeffs_value_x(2:end,:,1) - Fx(1:end-1,:) .* coeffs_value_x(1:end-1,:,2) ...
        - Dx(2:end,:) .* coeffs_derivative_x(2:end,:,1) + Dx(1:end-1,:) .* coeffs_derivative_x(1:end-1,:,2) ...
        + Fr(:,2:end) .* coeffs_value_r(:,2:end,1) - Fr(:,1:end-1) .* coeffs_value_r(:,1:end-1,2) ...
        - Dr(:,2:end) .* coeffs_derivative_r(:,2:end,1) + Dr(:,1:end-1) .* coeffs_derivative_r(:,1:end-1,2) ...
        - src_lin .* volume;

    % aW
    coeff(:,:,2) = Fx(1:end-1,:) .* coeffs_value_x(1:end-1,:,1) - Dx(1:end-1,:) .* coeffs_derivative_x(1:end-1,:,1);

    % aE
    coeff(:,:,3) = -Fx(2:end,:) .* coeffs_value_x(2:end,:,2) + Dx(2:end,:) .* coeffs_derivative_x(2:end,:,2);

    % aS
    coeff(:,:,4) = Fr(:,1:end-1) .* coeffs_value_r(:,1:end-1,1) - Dr(:,1:end-1) .* coeffs_derivative_r(:,1:end-1,1);

    % aN
    coeff(:,:,5) = -Fr(:,2:end) .* coeffs_value_r(:,2:end,2) + Dr(:,2:end) .* coeffs_derivative_r(:,2:end,2);

    % b
    coeff(:,:,6) = Fx(1:end-1,:) .* coeffs_value_x(1:end-1,:,3) - Fx(2:end,:) .* coeffs_value_x(2:end,:,3) ...
        - Dx(1:end-1,:) .* coeffs_derivative_x(1:end-1,:,3) + Dx(2:end,:) .* coeffs_derivative_x(2:end,:,3) ...
        + Fr(:,1:end-1) .* coeffs_value_r(:,1:end-1,3) - Fr(:,2:end) .* coeffs_value_r(:,2:end,3) ...
        - Dr(:,1:end-1) .* coeffs_derivative_r(:,1:end-1,3) + Dr(:,2:end) .* coeffs_derivative_r(:,2:end,3) ...
        + src .* volume;
end