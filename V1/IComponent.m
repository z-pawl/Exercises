classdef (Abstract) IComponent < handle & matlab.mixin.Heterogeneous
    methods (Abstract)
        update_properties(obj)
        
        converged = iterate(obj, coeff)

        res = calculate_residual(obj, coeff)
        converged = has_converged(obj)
        delete_residual_history(obj)
    end
end