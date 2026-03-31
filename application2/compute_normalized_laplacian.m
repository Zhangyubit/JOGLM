function Ln = compute_normalized_laplacian(bptG, norm_type)
% compute_normalized_laplacian - Compute normalized graph Laplacian matrix
%
% This function computes either the symmetric or asymmetric normalized
% Laplacian matrix for a given graph (or oversampled bipartite graph).
%
% Inputs:
%   bptG      - Adjacency matrix of the graph. It can be either:
%               (i) a matrix, or
%               (ii) a cell array where the adjacency matrix is stored in bptG{1}
%
%   norm_type - Type of normalization:
%               'sym'  : symmetric normalized Laplacian
%               'asym' : asymmetric (random-walk) normalized Laplacian
%
% Output:
%   Ln        - Normalized graph Laplacian matrix
%
% Notes:
%   - For the symmetric case:    Ln = I - D^{-1/2} A D^{-1/2}
%   - For the asymmetric case:  Ln = I - D^{-1} A
%   - Zero-degree nodes are handled by replacing zero entries with one
%     to avoid division-by-zero issues.

    % Extract adjacency matrix if input is a cell array
    if isa(bptG, 'cell')
        A = bptG{1};
    else
        A = bptG;
    end

    switch norm_type
        case 'sym'  % Symmetric normalized Laplacian
            d = sum(A, 2);
            d(d == 0) = 1;  % Avoid division by zero
            
            D_inv_sqrt = diag(d .^ (-0.5));
            An = D_inv_sqrt * A * D_inv_sqrt;
            
            % Enforce symmetry for numerical stability
            An = 0.5 * (An + An');
            
            Ln = eye(length(A)) - An;

        case 'asym' % Asymmetric (random-walk) normalized Laplacian
            % Enforce symmetry of adjacency (optional but improves stability)
            A = 0.5 * (A + A');
            
            d = sum(A, 2);
            d(d == 0) = 1;  % Avoid division by zero
            
            D_inv = diag(d .^ (-1));
            An = D_inv * A;
            
            Ln = eye(length(A)) - An;

        otherwise
            error('Unknown normalization type: norm_type must be "sym" or "asym".');
    end
end
