classdef Simulation2D < handle
    properties
        components (:,1) IComponent
        component_indexes (1,1) dictionary

        noi (1,1) int32
        max_iter (1,1) int32

        grid Grid2D {mustBeScalarOrEmpty}
    end
    methods
        function obj = update_properties(obj)
        end
        function ok = check_dependencies(obj)
            ok = all(arrayfun(@(c) c.check_dependencies(), obj.components));
        end
        function solve(obj)
            converged = zeros(size(obj.components));
            while ~all(converged) && obj.noi<=obj.max_iter
                obj = obj.update_properties();
                converged=arrayfun(@(c) c.iterate(), obj.components);
            end
        end
    end
end