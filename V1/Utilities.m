classdef Utilities
    methods (Static)
        % Linearly interpolates values from cell centers to cell faces along a specified direction
        function interpolated = interpolate_center2faces(val_at_center, cell_sizes, direction)
            % Large constant used to ensure that the interpolated value at
            % the boundary is calculated using only the value within the
            % domain
            INF = 1e25;

            % Changing the size of the cell sizes to be the same as the
            % val_at_center
            cell_sizes = paddata(cell_sizes,size(val_at_center),"Pattern","edge");
            
            % Padding size
            sz = size(val_at_center);
            sz(direction)=sz(direction)+1;

            % Padded values and cell sizes
            val_leading = paddata(val_at_center,sz,"Side","leading");
            val_trailing = paddata(val_at_center,sz,"Side","trailing");
            size_leading = paddata(cell_sizes,sz,"Side","leading","FillValue",INF);
            size_trailing = paddata(cell_sizes,sz,"Side","trailing","FillValue",INF);

            % Linear interpolation
            interpolated = (val_leading.*size_trailing + val_trailing.*size_leading) ...
                ./ (size_leading + size_trailing);
        end
        function s = sum_cell(cell_arr)
            s = 0;
            for i = 1:numel(cell_arr)
                s = s + cell_arr{i};
            end
        end
    end
end