function local_residual = calculate_local_residual(coeff, field)
    sz = size(field);

    % A column and a row of zeros used to pad the data in order to make sure that there is no index out of bounds error
    col0 = zeros(sz(1),1);
    row0 = zeros(1,sz(2));

    % r = aW*TW+aE*TE+aS*TS+aN*TN+b-aP*TP
    local_residual = [row0; coeff(2:end,:,1) .* field(1:end-1,:)] + [coeff(1:end-1,:,2) .* field(2:end,:); row0] ... 
        + [col0 coeff(:,2:end,3) .* field(:,1:end-1)] + [coeff(:,1:end-1,4) .* field(:,2:end) col0] ...
        + coeff(:,:,6) - coeff(:,:,5) .* field;
end