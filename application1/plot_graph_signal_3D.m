function plot_graph_signal_3D(G, f)
% plot_graph_signal_3D - 3D visualization of graph signals
%
% This function visualizes a graph signal in a 3D representation where:
%   - The graph structure is shown as edges on the 2D plane.
%   - The signal amplitude is represented as vertical bars at each node.
%   - Node colors encode signal intensity using a colormap.
%
% Inputs:
%   G - Graph structure containing:
%         G.coords : node coordinates (N x 2)
%         G.W      : adjacency matrix
%   f - Graph signal vector of size N x 1
%
% Description:
%   - Edges are plotted as thin gray lines on the base plane (z = 0).
%   - Each node is represented by a vertical bar whose height equals f(i).
%   - Node colors provide an additional encoding of signal magnitude.
%
% Requirements:
%   - G must contain node coordinates (G.coords)
%   - length(f) must match the number of nodes G.N

    %% Input validation
    if ~isfield(G, 'coords')
        error('Graph G must contain node coordinates (G.coords).');
    end
    if length(f) ~= G.N
        error('Signal length must match number of graph nodes.');
    end

    %% Extract node coordinates
    xpos = G.coords(:, 1);
    ypos = G.coords(:, 2);

    %% Initialize figure
    figure('Units', 'pixels', 'Position', [100, 100, 1024, 768]); 
    hold on;

    %% ------------------------------------------------------------
    %% Plot graph edges (base layer)
    %% ------------------------------------------------------------
    % Use upper triangular part to avoid duplicate edges
    [row, col] = find(triu(G.W));
    
    for i = 1:length(row)
        plot3([xpos(row(i)), xpos(col(i))], ...
              [ypos(row(i)), ypos(col(i))], ...
              [0, 0], ...
              'Color', [0.6, 0.6, 0.6], 'LineWidth', 0.5);
    end

    %% ------------------------------------------------------------
    %% Plot signal as vertical bars
    %% ------------------------------------------------------------
    for i = 1:G.N
        plot3([xpos(i), xpos(i)], ...
              [ypos(i), ypos(i)], ...
              [0, f(i)], ...
              'Color', [0.3, 0.8, 0.3], 'LineWidth', 3); % light green
    end

    %% ------------------------------------------------------------
    %% Plot nodes with color-coded signal values
    %% ------------------------------------------------------------
    scatter3(xpos, ypos, zeros(G.N,1), 70, f, 'filled');

    %% Custom colormap (blue → red gradient)
    nColors = 64;
    cmap = [linspace(0,1,nColors)', ...   % Red channel
            linspace(0,0,nColors)', ...   % Green channel
            linspace(1,0,nColors)'];     % Blue channel

    colormap(cmap);
    colorbar;   % Display color scale

    %% ------------------------------------------------------------
    %% Visualization settings
    %% ------------------------------------------------------------
    zlabel('Signal Value','FontSize',30,'FontWeight','bold');
    
    view(45, 37); 
    grid on; 
    axis tight; 
    box on;

    % Axis formatting
    set(gca,'FontSize',28,'FontWeight','bold'); 
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide x/y axes
    
    xlabel(''); 
    ylabel('');

    % Optional: hide colorbar if needed
    colorbar off;

end