function [GFT, H0, H1, G0, G1, connect_bpt, Colorednodes] = ogfb_filters(G, filterlen, norm_type)
% ogfb_filters - Construct oversampled graph filter banks (OGFB)
%
% This function designs a two-channel graph filter bank on an oversampled
% bipartite graph obtained via graph lifting and coloring-based partitioning.
%
% Inputs:
%   G           - Graph structure (GSPBox format, e.g., gsp_ring, gsp_erdos_renyi)
%   filterlen   - Filter length (recommended to be even)
%   norm_type   - Type of Laplacian normalization:
%                 'sym'  : symmetric normalized Laplacian
%                 'asym' : random-walk normalized Laplacian
%
% Outputs:
%   GFT          - Graph Fourier transform matrix (inverse eigenvector matrix)
%   H0, H1       - Analysis filters (low-pass and high-pass)
%   G0, G1       - Synthesis filters (low-pass and high-pass)
%   connect_bpt  - Indices of duplicated nodes (oversampling connections)
%   Colorednodes - Node partitions obtained from graph coloring
%
% Description:
%   1) The input graph is first colored using the DSATUR algorithm.
%   2) A bipartite decomposition is obtained via Harary decomposition.
%   3) An oversampled graph Laplacian model (OSGLM) is constructed by
%      duplicating selected nodes according to the coloring structure.
%   4) A normalized Laplacian is computed on the lifted (oversampled) graph.
%   5) Biorthogonal spectral kernels are used to construct analysis and
%      synthesis filter banks in the graph spectral domain.
%
% Key idea:
%   Oversampling is achieved via graph lifting:
%       G  →  G_os
%   where G_os contains duplicated nodes indexed by connect_bpt, enabling
%   redundant multiresolution representations.

    A = G.W;

    %% Step 1: Graph coloring

    % Apply DSATUR coloring algorithm
    F = DSATUR(A);

    % Obtain color partitions via Harary decomposition
    [~, ~, ~, Colorednodes] = harary_decomp(A, F);

    % Recover coloring index vector
    F = recover_colorednodes(Colorednodes);

    %% Step 2: Construct oversampled bipartite graph (OSGLM)

    % Create oversampled graph Laplacian model
    [bptG, Colorednodes, connect_bpt] = ...
        create_OSGLM(A, Colorednodes, 2^ceil(log2(max(F))));

    %% Step 3: Compute normalized Laplacian

    Ln = compute_normalized_laplacian(bptG{1}, norm_type);

    %% Step 4: Design spectral filter kernels

    % Split filter length
    N_lo = floor(filterlen / 2) + mod(filterlen, 2);
    N_hi = N_lo;

    % Generate biorthogonal kernels
    [lo_d, hi_d] = biorth_kernel(N_lo, N_hi);

    % Analysis filters
    h0 = @(x) polyval(lo_d, x);        % Low-pass
    h1 = @(x) polyval(hi_d, x);        % High-pass

    % Synthesis filters (spectral flipping)
    g0 = @(x) polyval(hi_d, 2 - x);
    g1 = @(x) polyval(lo_d, 2 - x);

    %% Step 5: Graph Fourier Transform (GFT)

    [GFT, V, Lambda] = gft(Ln);

    %% Step 6: Construct filtering operators

    % Analysis operators
    H0 = V * diag(h0(diag(Lambda))) * GFT;
    H1 = V * diag(h1(diag(Lambda))) * GFT;

    % Synthesis operators
    G0 = V * diag(g0(diag(Lambda))) * GFT;
    G1 = V * diag(g1(diag(Lambda))) * GFT;

end
