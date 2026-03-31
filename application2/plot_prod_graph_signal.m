function plot_prod_graph_signal(G_base, N_time, x)
% plot_prod_graph_signal - Visualization of signals on a product graph
%
% This function visualizes a graph signal defined on the product graph:
%       G = G_base ⊗ G_ring(N_time)
%
% The product graph is embedded in 3D space, where:
%   - The base graph (G_base) defines spatial connections (x-y plane)
%   - The time dimension is represented along the z-axis
%   - The signal is encoded using node colors
%
% Inputs:
%   G_base - Base graph structure containing:
%              G_base.W      : adjacency matrix
%              G_base.coords : node coordinates (N_base x 2)
%   N_time - Number of time samples (e.g., 4)
%   x      - Signal vector of length N_base * N_time
%
% Description:
%   - Each node in the product graph corresponds to (vertex, time).
%   - Spatial edges connect nodes within the same time layer.
%   - Temporal edges connect the same vertex across time (ring structure).
%   - The signal is visualized using color intensity at each node.

    %% Extract base graph information
    N_base     = G_base.N;
    coords_base = G_base.coords;
    W_base      = G_base.W;

    %% ============================================================
    %% Step 1: Construct product graph coordinates (3D embedding)
    %% ============================================================

    % Allocate coordinate matrix
    coords = zeros(N_base * N_time, 3);

    for i = 1:N_base
        for t = 1:N_time
            idx = (i-1)*N_time + t;

            % Spatial coordinates (x, y)
            coords(idx,1:2) = coords_base(i,:);

            % Temporal coordinate (z-axis)
            coords(idx,3) = t;
        end
    end

    %% ============================================================
    %% Step 2: Initialize figure
    %% ============================================================

    figure('Color','w','Units','pixels','Position',[100,100,1024,768]);
    hold on;

    %% ============================================================
    %% Step 3: Plot spatial edges (within each time layer)
    %% ============================================================

    for i = 1:N_base
        for j = i+1:N_base
            if W_base(i,j) > 0
                for t = 1:N_time
                    idx1 = (i-1)*N_time + t;
                    idx2 = (j-1)*N_time + t;

                    plot3([coords(idx1,1), coords(idx2,1)], ...
                          [coords(idx1,2), coords(idx2,2)], ...
                          [coords(idx1,3), coords(idx2,3)], ...
                          'Color', [0.6 0.6 0.6], 'LineWidth', 1.0);
                end
            end
        end
    end

    %% ============================================================
    %% Step 4: Plot temporal edges (ring structure)
    %% ============================================================

    for i = 1:N_base
        for t = 1:N_time
            idx1 = (i-1)*N_time + t;

            % Ring connection in time dimension
            idx2 = (i-1)*N_time + mod(t, N_time) + 1;

            plot3([coords(idx1,1), coords(idx2,1)], ...
                  [coords(idx1,2), coords(idx2,2)], ...
                  [coords(idx1,3), coords(idx2,3)], ...
                  'Color', [0.6 0.6 0.6], 'LineWidth', 1.0);
        end
    end

    %% ============================================================
    %% Step 5: Visualize signal using node colors
    %% ============================================================

    x_real = real(x(:));

    % Normalize signal for color mapping
    c_min = min(x_real);
    c_max = max(x_real);
    colors = (x_real - c_min) / (c_max - c_min + eps);

    scatter3(coords(:,1), coords(:,2), coords(:,3), ...
             100, colors, 'filled');

    colormap(parula);

    cb = colorbar;
    set(cb, 'FontSize', 16, 'FontWeight', 'bold');

    %% ============================================================
    %% Step 6: Visualization settings
    %% ============================================================

    view(45,30);   % Adjust viewing angle
    axis off;
    grid on;

    set(gca,'Color','none');
    set(gcf,'Color','w');

end