function F_new = recover_colorednodes(Colorednodes)
% recover_colorednodes - Recover node color labels from partitioned node sets
%
% This function reconstructs a color label vector from a cell array of
% node partitions obtained via graph coloring.
%
% Input:
%   Colorednodes - Cell array where each cell contains the indices of nodes
%                  belonging to a specific color class
%
% Output:
%   F_new        - Column vector of size N x 1, where F_new(i) denotes
%                  the color label assigned to node i
%
% Description:
%   - Each cell in Colorednodes corresponds to a color group.
%   - The function assigns the corresponding color index to each node.
%   - Empty cells are safely ignored to ensure robustness.

    %% Determine the maximum node index (graph size)

    % Identify non-empty cells
    non_empty = ~cellfun(@isempty, Colorednodes);

    if any(non_empty)
        % Compute maximum node index across all non-empty cells
        max_vals = cellfun(@(x) max(x(:)), ...
                           Colorednodes(non_empty), ...
                           'UniformOutput', true);
        N = max(max_vals);
    else
        % Handle case where all cells are empty
        N = 0;
    end

    %% Initialize color label vector
    F_new = zeros(N, 1);

    %% Assign color labels to nodes
    for color = 1:length(Colorednodes)
        nodes = Colorednodes{color};
        
        if ~isempty(nodes)
            F_new(nodes) = color;
        end
    end

end
