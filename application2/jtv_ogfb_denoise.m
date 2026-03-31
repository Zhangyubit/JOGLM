function X_rec = jtv_ogfb_denoise(G, X, filterlen, norm_type, sigma, thresh)
% jtv_ogfb_denoise - Joint time-vertex denoising using oversampled graph filter banks
%
% This function performs denoising of a time-vertex graph signal using an
% oversampled joint time-vertex graph filter bank (OGFB). Redundancy is
% introduced via graph lifting, and noise is suppressed by thresholding
% high-frequency components in both time and vertex domains.
%
% Input:
%   G          - Graph structure (from GSPBox)
%   X          - Graph signal matrix of size N × T
%   filterlen  - Length of the graph filters
%   norm_type  - Laplacian normalization type ('sym' or 'asym')
%   sigma      - Noise standard deviation
%   thresh     - Threshold scaling factor
%
% Output:
%   X_rec      - Reconstructed (denoised) graph signal

disp('***** Joint Time-Vertex Oversampled Graph Filter Bank Denoising *****')

%% Basic parameters
N = size(X,1);         % Number of vertices
T = size(X,2);         % Number of time samples

% Construct time graph (ring structure)
GT = gsp_ring(T);

%% Construct oversampled filter banks
[~, H0_T, H1_T, G0_T, G1_T, connect_bpt_T, Colorednodes_T] = ...
    ogfb_filters(GT, filterlen, norm_type);

[~, H0_G, H1_G, G0_G, G1_G, connect_bpt_G, Colorednodes_G] = ...
    ogfb_filters(G, filterlen, norm_type);


%% ================= OVERSAMPLED SIGNAL CONSTRUCTION =================
% Extend signal along time and vertex domains (graph lifting)

X_T_os = [X, X(:, connect_bpt_T)];          % Time-domain oversampling
X_os   = [X_T_os; X_T_os(connect_bpt_G, :)];% Vertex-domain oversampling

[N_os, T_os] = size(X_os);                  % Size of oversampled signal


%% ================= TIME DOMAIN PROCESSING =================

% Identify duplicated indices in bipartite structure
[~, OS1_T, ~] = intersect(connect_bpt_T, Colorednodes_T{5});
[~, OS2_T, ~] = intersect(connect_bpt_T, Colorednodes_T{6});

%% ----- Decomposition -----
X_wL_T = zeros(N_os, T_os); 
X_wH_T = zeros(N_os, T_os);

tempL_T = X_os * conj(H0_T); 
tempH_T = X_os * conj(H1_T);

% Sparse assignment including duplicated nodes
X_wL_T(:, Colorednodes_T{5}) = tempL_T(:, Colorednodes_T{5});
X_wL_T(:, T + OS2_T)         = tempL_T(:, T + OS2_T);

X_wH_T(:, Colorednodes_T{6}) = tempH_T(:, Colorednodes_T{6});
X_wH_T(:, T + OS1_T)         = tempH_T(:, T + OS1_T);

%% ----- Thresholding (high-frequency components) -----
thresh_T = thresh * sigma;

X_wH_T = X_wH_T .* (abs(X_wH_T) > thresh_T);

%% ----- Reconstruction -----
tempX_hatL_T = X_wL_T * conj(G0_T);
tempX_hatH_T = X_wH_T * conj(G1_T);

X_hat_T = tempX_hatL_T + tempX_hatH_T;


%% ================= VERTEX DOMAIN PROCESSING =================

[~, OS1_G, ~] = intersect(connect_bpt_G, Colorednodes_G{5});
[~, OS2_G, ~] = intersect(connect_bpt_G, Colorednodes_G{6});

%% ----- Decomposition -----
X_wL_G = zeros(N_os, T_os); 
X_wH_G = zeros(N_os, T_os);

tempL_G = H0_G * X_hat_T; 
tempH_G = H1_G * X_hat_T;

% Sparse assignment including duplicated nodes
X_wL_G(Colorednodes_G{5}, :) = tempL_G(Colorednodes_G{5}, :);
X_wL_G(N + OS2_G, :)         = tempL_G(N + OS2_G, :);

X_wH_G(Colorednodes_G{6}, :) = tempH_G(Colorednodes_G{6}, :);
X_wH_G(N + OS1_G, :)         = tempH_G(N + OS1_G, :);

%% ----- Thresholding (high-frequency components) -----
thresh_G = thresh * sigma;

X_wH_G = X_wH_G .* (abs(X_wH_G) > thresh_G);

%% ----- Reconstruction -----
tempX_hatL_G = G0_G * X_wL_G;
tempX_hatH_G = G1_G * X_wH_G;

X_hat_G = tempX_hatL_G + tempX_hatH_G;


%% ================= DOWNSAMPLING (FOLDING BACK) =================
% Merge duplicated nodes and recover original signal size

% Vertex domain folding
X_hat_G_final = X_hat_G(1:N, :);
X_hat_G_final(connect_bpt_G, :) = ...
    0.5 * (X_hat_G(connect_bpt_G, :) + X_hat_G(N+1:end, :));

% Time domain folding
X_hat_final = X_hat_G_final(:, 1:T);
X_hat_final(:, connect_bpt_T) = ...
    0.5 * (X_hat_G_final(:, connect_bpt_T) + X_hat_G_final(:, T+1:end));

%% Final reconstructed signal
X_rec = X_hat_final;

end