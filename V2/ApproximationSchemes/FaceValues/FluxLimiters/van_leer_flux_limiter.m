function flux_limiter = van_leer_flux_limiter(derivative_face, derivative_upwind_face)
    arguments
        derivative_face double
        derivative_upwind_face double 
    end

    if derivative_face ~= 0
        r = derivative_upwind_face / derivative_face;
        if r >= 0
            flux_limiter = 2 * r / (1 + r);
        else
            flux_limiter = 0;
        end
    else
        flux_limiter = 0;
    end
end