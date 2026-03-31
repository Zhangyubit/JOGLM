function plot_oversampled_bipartite_graph(G)
% plot_oversampled_bipartite_graph - Visualization of oversampled bipartite graph
%
% This function visualizes the structure of an oversampled bipartite graph
% obtained via graph lifting and coloring-based partitioning.
%
% Inputs:
%   G - Graph structure containing:
%         G.W      : adjacency matrix
%         G.coords : node coordinates (2D)
%
% Description:
%   - The original graph is first decomposed using graph coloring and
%     Harary decomposition.
%   - An oversampled bipartite graph is constructed via node duplication.
%   - The resulting graph is visualized in 3D:
%         • Bottom layer  : original nodes
%         • Top layer     : duplicated nodes (oversampling)
%   - Nodes are classified into:
%         • Low-pass nodes  (blue squares)
%         • High-pass nodes (red circles)
%   - Edges include:
%         • Intra-layer edges (original graph structure)
%         • Inter-layer edges (oversampling connections)

    A = G.W;
    coords = G.coords;
    N = size(A,1);

    %% ============================================================
    %% Step 1: Graph coloring and Harary decomposition
    %% ============================================================

    % Apply DSATUR graph coloring
    F = DSATUR(A);

    % Obtain color partitions
    [~, ~, ~, Colorednodes] = harary_decomp(A, F);

    % Recover color index vector
    F = recover_colorednodes(Colorednodes);

    %% ============================================================
    %% Step 2: Construct oversampled bipartite graph (OSGLM)
    %% ============================================================

    % Generate oversampled graph via node duplication
    [bptG, Colorednodes_bpt, connect_bpt] = ...
        create_OSGLM(A, Colorednodes, 2^ceil(log2(max(F))));

    A_os = bptG{1};   % Adjacency matrix of oversampled graph

    %% ============================================================
    %% Step 3: Extend node coordinates (3D lifting)
    %% ============================================================

    dz = 0.0001;   % Height offset between layers

    % Construct 3D coordinates:
    %   First N nodes: original layer (z = 0)
    %   Next  N nodes: duplicated layer (z = dz)
    coords_os = [coords, zeros(N,1); ...
                 coords, dz*ones(N,1)];

    %% ============================================================
    %% Step 4: Identify low-pass and high-pass nodes
    %% ============================================================

    % Identify duplicated node indices in each partition
    [~, OS1] = intersect(connect_bpt, Colorednodes_bpt{5});
    [~, OS2] = intersect(connect_bpt, Colorednodes_bpt{6});

    % Original layer
    lowpass1  = Colorednodes_bpt{5};
    highpass1 = Colorednodes_bpt{6};

    % Duplicated layer
    lowpass2  = N + OS2;
    highpass2 = N + OS1;

    % Combine node sets
    lowpass_nodes  = [lowpass1;  lowpass2];
    highpass_nodes = [highpass1; highpass2];

    %% ============================================================
    %% Step 5: Extract edges
    %% ============================================================

    [edge_i, edge_j] = find(A_os > 0);

    edges = [edge_i, edge_j];
    edges = edges(edge_i < edge_j, :);   % Remove duplicate edges

    %% ============================================================
    %% Step 6: Visualization
    %% ============================================================

    figure('Color', 'w', 'Position', [100, 100, 1024, 768]);
    hold on; axis off; view(3); grid on;

    xlabel('X'); ylabel('Y'); zlabel('Z');

    %% ---- Plot edges ----
    for k = 1:size(edges,1)
        i = edges(k,1); 
        j = edges(k,2);

        xi = coords_os(i,1); yi = coords_os(i,2); zi = coords_os(i,3);
        xj = coords_os(j,1); yj = coords_os(j,2); zj = coords_os(j,3);

        % Check if edge is inter-layer (vertical connection)
        if abs(zi - zj) > 0
            % Inter-layer edges (oversampling links)
            plot3([xi, xj], [yi, yj], [zi, zj], ...
                  'Color', [0, 1, 0, 0.2], ... % light transparent green
                  'LineWidth', 1.2);
        else
            % Intra-layer edges (original graph structure)
            plot3([xi, xj], [yi, yj], [zi, zj], ...
                  'Color', [0.5, 0.5, 0.5], ...
                  'LineWidth', 0.5);
        end
    end

    %% ---- Plot nodes ----

    % Low-pass nodes: blue squares
    plot3(coords_os(lowpass_nodes,1), ...
          coords_os(lowpass_nodes,2), ...
          coords_os(lowpass_nodes,3), ...
          's', 'MarkerEdgeColor', 'b', ...
          'MarkerFaceColor', 'none', ...
          'MarkerSize', 8, 'LineWidth', 1.5);

    % High-pass nodes: red circles
    plot3(coords_os(highpass_nodes,1), ...
          coords_os(highpass_nodes,2), ...
          coords_os(highpass_nodes,3), ...
          'o', 'MarkerEdgeColor', 'r', ...
          'MarkerFaceColor', 'none', ...
          'MarkerSize', 8, 'LineWidth', 1.5);

    %% Final view
    view([20, 70]);

end