function [GFT, H0, H1, G0, G1, Colorednodes] = gfb_filters(G, filterlen, norm_type)
% gfb_filters - Construct graph filter bank on an oversampled bipartite graph
%
% This function designs a two-channel graph filter bank by first constructing
% a bipartite subgraph via graph coloring, and then defining spectral filters
% based on polynomial kernels in the graph spectral domain.
%
% Inputs:
%   G          - Graph structure from GSPBox (e.g., gsp_ring, gsp_erdos_renyi)
%   filterlen  - Filter length (recommended to be even)
%   norm_type  - Type of Laplacian normalization:
%                'sym'  : symmetric normalized Laplacian
%                'asym' : random-walk (asymmetric) normalized Laplacian
%
% Outputs:
%   GFT          - Graph Fourier transform matrix (i.e., V^{-1})
%   H0, H1       - Analysis filters (low-pass and high-pass)
%   G0, G1       - Synthesis filters (low-pass and high-pass)
%   Colorednodes - Cell array of node partitions obtained via graph coloring,
%                  used for bipartite construction and sparse processing
%
% Description:
%   1) The graph is first colored using the DSATUR algorithm.
%   2) Nodes are partitioned into two disjoint sets by merging odd- and even-indexed color classes.
%   3) A bipartite subgraph is constructed from these two node sets.
%   4) The normalized Laplacian is computed on the bipartite graph.
%   5) Polynomial filter kernels are applied in the spectral domain to construct
%      analysis and synthesis filter banks.
%
% Notes:
%   - The bipartite structure is essential for enabling perfect reconstruction.
%   - Filters are designed via biorthogonal polynomial kernels.
%   - Colored node groups are preserved for later sparse sampling operations.

    A = G.W;

    %% Graph coloring and bipartite partition construction
    
    % Apply DSATUR graph coloring algorithm
    F = DSATUR(A);
    
    % Perform Harary decomposition to obtain color classes
    [~,~,~, Colorednodes] = harary_decomp(A, F);

    % Initialize two merged node groups (for bipartite graph)
    Colorednodes{5} = [];
    Colorednodes{6} = [];

    % Merge odd-indexed and even-indexed color classes
    for c = 1:numel(Colorednodes)
        if mod(c,2) == 1
            Colorednodes{5} = union(Colorednodes{5}, Colorednodes{c});
        else
            Colorednodes{6} = union(Colorednodes{6}, Colorednodes{c});
        end
    end

    %% Construct bipartite adjacency matrix
    
    A_bpt = zeros(size(A));
    
    % Keep only inter-set connections
    A_bpt(Colorednodes{5}, Colorednodes{6}) = A(Colorednodes{5}, Colorednodes{6});
    A_bpt(Colorednodes{6}, Colorednodes{5}) = A(Colorednodes{6}, Colorednodes{5});
    
    % Ensure symmetry
    A_bpt = (A_bpt + A_bpt.') / 2;

    %% Compute normalized Laplacian
    
    Ln = compute_normalized_laplacian(A_bpt, norm_type);

    %% Design spectral filter kernels
    
    % Split filter length into low-pass and high-pass components
    N_lo = floor(filterlen / 2) + mod(filterlen, 2);
    N_hi = N_lo;

    % Generate biorthogonal polynomial kernels
    [lo_d, hi_d] = biorth_kernel(N_lo, N_hi);

    % Define analysis filters (spectral domain)
    h0 = @(x) polyval(lo_d, x);          % Low-pass analysis filter
    h1 = @(x) polyval(hi_d, x);          % High-pass analysis filter

    % Define synthesis filters (spectral domain)
    g0 = @(x) polyval(hi_d, 2 - x);      % Low-pass synthesis filter
    g1 = @(x) polyval(lo_d, 2 - x);      % High-pass synthesis filter

    %% Graph Fourier Transform (GFT)
    
    [GFT, V, Lambda] = gft(Ln);  % V: eigenvectors, Lambda: eigenvalues

    %% Construct filtering operators
    
    % Analysis operators
    H0 = V * diag(h0(diag(Lambda))) * GFT;
    H1 = V * diag(h1(diag(Lambda))) * GFT;

    % Synthesis operators
    G0 = V * diag(g0(diag(Lambda))) * GFT;
    G1 = V * diag(g1(diag(Lambda))) * GFT;

end