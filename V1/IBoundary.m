classdef (Abstract) IBoundary < handle & matlab.mixin.Heterogeneous
    properties (Abstract)
        name (1,1) string
        field_sz double
        boundaries (2,:) 
    end
    methods (Abstract)
        set_boundary_condtion(obj, boundary_index, boundary_side)
        coeff = apply_boundary_conditions(obj, coeff)
    end
end