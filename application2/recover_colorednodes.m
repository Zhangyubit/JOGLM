function F_new = recover_colorednodes(Colorednodes)
% recover_colorednodes - Convert color partitions to a node label vector
%
% This function reconstructs a node-wise color label vector from a cell
% array of node partitions obtained via graph coloring.
%
% Input:
%   Colorednodes - Cell array where each cell contains indices of nodes
%                  belonging to a specific color class
%
% Output:
%   F_new        - N x 1 vector, where F_new(i) is the color label of node i
%
% Description:
%   - Each cell corresponds to one color class.
%   - The function assigns the corresponding color index to each node.
%   - Empty cells are ignored to ensure robustness.

    %% Determine graph size (maximum node index)
    non_empty = ~cellfun(@isempty, Colorednodes);

    if any(non_empty)
        max_vals = cellfun(@(x) max(x(:)), ...
                           Colorednodes(non_empty), ...
                           'UniformOutput', true);
        N = max(max_vals);
    else
        N = 0;   % All cells are empty
    end

    %% Initialize label vector
    F_new = zeros(N, 1);

    %% Assign color labels
    for color = 1:length(Colorednodes)
        nodes = Colorednodes{color};
        if ~isempty(nodes)
            F_new(nodes) = color;
        end
    end

end