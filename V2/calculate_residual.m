function residual = calculate_residual(coeff, field)
    % r = aW*TW+aE*TE+aS*TS+aN*TN+b-aP*TP
    local_residual = calculate_local_residual(coeff, field);
    
    % s_f = || [aP1*TP1; aP2*TP2; aP3*TP3; .... aPN-1*TPN-1; aPN*TPN] ||
    scaling_factor = norm(coeff(:,:,5) .* field, 1);
    
    % residual = || r || / s_f
    residual = norm(local_residual, 1) / scaling_factor;
end