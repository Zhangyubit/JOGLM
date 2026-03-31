function X_rec = jtv_gfb_denoise(G, X, filterlen, norm_type, sigma, thresh)
% jtv_gfb_denoise - Joint time-vertex denoising using critically sampled graph filter banks
%
% This function performs denoising of a time-vertex graph signal using a 
% critically sampled joint time-vertex graph filter bank (GFB). The method
% consists of two stages: time-domain processing followed by vertex-domain
% processing. Thresholding is applied to suppress high-frequency noise.
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

disp('***** Joint Time-Vertex Graph Filter Bank Denoising *****')

%% Basic parameters
N = size(X,1);         % Number of vertices
T = size(X,2);         % Number of time samples

% Construct time graph (ring topology)
GT = gsp_ring(T);

%% Construct filter banks for time and vertex domains
[~, H0_T, H1_T, G0_T, G1_T, Colorednodes_T] = gfb_filters(GT, filterlen, norm_type);
[~, H0_G, H1_G, G0_G, G1_G, Colorednodes_G] = gfb_filters(G,  filterlen, norm_type);


%% ================= TIME DOMAIN PROCESSING =================
% Decomposition → Thresholding → Reconstruction

%% ----- Decomposition -----
X_wL_T = zeros(N, T); 
X_wH_T = zeros(N, T);

tempL_T = X * conj(H0_T); 
tempH_T = X * conj(H1_T);

% Sparse assignment based on bipartite partition
X_wL_T(:, Colorednodes_T{5}) = tempL_T(:, Colorednodes_T{5});
X_wH_T(:, Colorednodes_T{6}) = tempH_T(:, Colorednodes_T{6});

%% ----- Thresholding (high-frequency components) -----
thresh_T = thresh * sigma;

% Only threshold high-frequency coefficients
X_wH_T = X_wH_T .* (abs(X_wH_T) > thresh_T);

%% ----- Reconstruction -----
tempX_hatL_T = X_wL_T * conj(G0_T);
tempX_hatH_T = X_wH_T * conj(G1_T);

X_hat_T = tempX_hatL_T + tempX_hatH_T;


%% ================= VERTEX DOMAIN PROCESSING =================
% Decomposition → Thresholding → Reconstruction

%% ----- Decomposition -----
X_wL_G = zeros(N, T); 
X_wH_G = zeros(N, T);

tempL_G = H0_G * X_hat_T; 
tempH_G = H1_G * X_hat_T;

% Sparse assignment based on bipartite partition
X_wL_G(Colorednodes_G{5}, :) = tempL_G(Colorednodes_G{5},:);
X_wH_G(Colorednodes_G{6}, :) = tempH_G(Colorednodes_G{6},:);

%% ----- Thresholding (high-frequency components) -----
thresh_G = thresh * sigma;

% Only threshold high-frequency coefficients
X_wH_G = X_wH_G .* (abs(X_wH_G) > thresh_G);

%% ----- Reconstruction -----
tempX_hatL_G = G0_G * X_wL_G;
tempX_hatH_G = G1_G * X_wH_G;

X_hat_G = tempX_hatL_G + tempX_hatH_G;

%% Final reconstructed signal
X_rec = X_hat_G;

end
