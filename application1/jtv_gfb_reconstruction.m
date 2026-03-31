function [hatX, X_rec, MSE, NMSE, SNR] = jtv_gfb_reconstruction(G, X)
% jtv_gfb_reconstruction - Joint time-vertex reconstruction using critically sampled graph filter banks
%
% This function performs joint time-vertex graph signal reconstruction using
% a separable two-stage graph filter bank framework:
%   1) Temporal-domain filtering (time graph)
%   2) Vertex-domain filtering (graph structure)
%
% The method is based on critically sampled graph filter banks constructed
% via bipartite graph decomposition.
%
% Inputs:
%   G     - Graph structure (GSPBox format)
%   X     - Graph signal matrix of size N x T (N nodes, T time instances)
%
% Outputs:
%   hatX  - Joint time-vertex spectral representation (oversampled JFT)
%   X_rec - Reconstructed graph signal (N x T)
%   MSE   - Mean Squared Error
%   NMSE  - Normalized Mean Squared Error
%   SNR   - Signal-to-Noise Ratio (dB)
%
% Description:
%   - The signal is first processed along the time dimension using a ring graph.
%   - Then, the intermediate result is processed along the vertex domain.
%   - Each stage includes analysis (decomposition), masking (sampling),
%     and synthesis (reconstruction).
%   - Sparse assignments are enforced using colored node partitions.

%% Parameter settings
filterlen = 12;
norm_type = 'sym';

T = size(X, 2);
N = G.N;

%% Construct temporal graph
GT = gsp_ring(T);

X_or = X;                      % Original signal
f_hat_or = vec(X_or);          % Vectorized reference signal

%% ============================================================
%% Step 0: Construct joint time-vertex filter banks
%% ============================================================

% Time-domain filter bank
[F_T, H0_T, H1_T, G0_T, G1_T, Colorednodes_T] = gfb_filters(GT, filterlen, norm_type);

% Vertex-domain filter bank
[F_G, H0_G, H1_G, G0_G, G1_G, Colorednodes_G] = gfb_filters(G, filterlen, norm_type);

%% ============================================================
%% Step 1: Joint Time-Vertex Fourier Transform (JFT)
%% ============================================================

% Compute joint spectral representation
hatX = F_G * X_or * conj(F_T);

%% ============================================================
%% Step 2: Temporal domain processing
%% ============================================================

% Initialize low-pass and high-pass channels
X_wL_T = zeros(N, T); 
X_wH_T = zeros(N, T);

% Analysis filtering (time direction)
tempL_T = X_or * conj(H0_T); 
tempH_T = X_or * conj(H1_T);

% Sparse assignment based on bipartite partition
X_wL_T(:, Colorednodes_T{5}) = tempL_T(:, Colorednodes_T{5});
X_wH_T(:, Colorednodes_T{6}) = tempH_T(:, Colorednodes_T{6});

% Synthesis (time reconstruction)
tempX_hatL_T = X_wL_T * conj(G0_T);
tempX_hatH_T = X_wH_T * conj(G1_T);

% Combine channels
X_hat_T = tempX_hatL_T + tempX_hatH_T;

%% ============================================================
%% Step 3: Vertex domain processing
%% ============================================================

% Initialize channels
X_wL_G = zeros(N, T); 
X_wH_G = zeros(N, T);

% Analysis filtering (vertex domain)
tempL_G = H0_G * X_hat_T; 
tempH_G = H1_G * X_hat_T;

% Sparse assignment based on graph bipartition
X_wL_G(Colorednodes_G{5}, :) = tempL_G(Colorednodes_G{5},:);
X_wH_G(Colorednodes_G{6}, :) = tempH_G(Colorednodes_G{6},:);

% Synthesis (vertex reconstruction)
tempX_hatL_G = G0_G * X_wL_G;
tempX_hatH_G = G1_G * X_wH_G;

% Combine channels
X_hat_G = tempX_hatL_G + tempX_hatH_G;

%% Final reconstructed signal
X_rec = X_hat_G;

%% ============================================================
%% Performance evaluation
%% ============================================================

f_hat = vec(X_rec);

MSE  = mean((f_hat_or - f_hat).^2);
NMSE = calculate_nmse(f_hat_or, f_hat);

P_rms = mean(f_hat_or.^2);
SNR   = 10 * log10(P_rms / MSE);

end