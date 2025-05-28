classdef (Abstract) IComponent < handle & matlab.mixin.Heterogeneous
    methods (Abstract)
        converged = iterate(obj)
    end
end