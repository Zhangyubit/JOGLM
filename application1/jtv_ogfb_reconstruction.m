function [hatX_os, X_rec, MSE, NMSE, SNR] = jtv_ogfb_reconstruction(G, X)
% jtv_ogfb_reconstruction - Joint time-vertex reconstruction using oversampled graph filter banks
%
% This function performs joint time-vertex graph signal reconstruction using
% an oversampled graph filter bank (OGFB) framework. The method introduces
% redundancy via bipartite expansion in both time and vertex domains,
% enabling improved reconstruction performance and flexible representations.
%
% Inputs:
%   G     - Graph structure (GSPBox format)
%   X     - Graph signal matrix of size N x T
%
% Outputs:
%   hatX_os - Joint spectral representation (oversampled JFT)
%   X_rec   - Reconstructed graph signal (N x T)
%   MSE     - Mean Squared Error
%   NMSE    - Normalized Mean Squared Error
%   SNR     - Signal-to-Noise Ratio (dB)
%
% Description:
%   - The signal is first lifted into an oversampled joint domain by duplicating
%     selected nodes and time indices based on bipartite connections.
%   - A separable filter bank is applied:
%         (1) Time-domain filtering
%         (2) Vertex-domain filtering
%   - Reconstruction is followed by a merging (downsampling) step to remove redundancy.
%
% Key idea:
%   Oversampling is achieved by augmenting the signal:
%       X_os = [ X        X(:, connect_T);
%                X(connect_G, :)  ...      ]
%   which enables redundant multiresolution representations.

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
%% Step 0: Construct oversampled graph filter banks
%% ============================================================

% Time-domain OGFB
[F_T, H0_T, H1_T, G0_T, G1_T, connect_bpt_T, Colorednodes_T] = ...
    ogfb_filters(GT, filterlen, norm_type);

% Vertex-domain OGFB
[F_G, H0_G, H1_G, G0_G, G1_G, connect_bpt_G, Colorednodes_G] = ...
    ogfb_filters(G, filterlen, norm_type);

%% ============================================================
%% Step 1: Oversampled joint signal construction
%% ============================================================

% Time-domain oversampling (column duplication)
X_T_os = [X, X(:, connect_bpt_T)];

% Vertex-domain oversampling (row duplication)
X_os = [X_T_os; X_T_os(connect_bpt_G, :)];

[N_os, T_os] = size(X_os);

%% ============================================================
%% Step 2: Oversampled Joint Fourier Transform (JFT)
%% ============================================================

hatX_os = F_G * X_os * conj(F_T);

%% ============================================================
%% Step 3: Temporal domain processing
%% ============================================================

% Identify duplicated indices within bipartite sets
[~, OS1_T, ~] = intersect(connect_bpt_T, Colorednodes_T{5});
[~, OS2_T, ~] = intersect(connect_bpt_T, Colorednodes_T{6});

% Initialize channels
X_wL_T = zeros(N_os, T_os); 
X_wH_T = zeros(N_os, T_os);

% Analysis filtering
tempL_T = X_os * conj(H0_T); 
tempH_T = X_os * conj(H1_T);

% Sparse assignment including duplicated components
X_wL_T(:, Colorednodes_T{5}) = tempL_T(:, Colorednodes_T{5});
X_wL_T(:, T + OS2_T)         = tempL_T(:, T + OS2_T);

X_wH_T(:, Colorednodes_T{6}) = tempH_T(:, Colorednodes_T{6});
X_wH_T(:, T + OS1_T)         = tempH_T(:, T + OS1_T);

% Synthesis
tempX_hatL_T = X_wL_T * conj(G0_T);
tempX_hatH_T = X_wH_T * conj(G1_T);

X_hat_T = tempX_hatL_T + tempX_hatH_T;

%% ============================================================
%% Step 4: Vertex domain processing
%% ============================================================

% Identify duplicated indices
[~, OS1_G, ~] = intersect(connect_bpt_G, Colorednodes_G{5});
[~, OS2_G, ~] = intersect(connect_bpt_G, Colorednodes_G{6});

% Initialize channels
X_wL_G = zeros(N_os, T_os); 
X_wH_G = zeros(N_os, T_os);

% Analysis filtering
tempL_G = H0_G * X_hat_T; 
tempH_G = H1_G * X_hat_T;

% Sparse assignment including duplicated nodes
X_wL_G(Colorednodes_G{5}, :) = tempL_G(Colorednodes_G{5}, :);
X_wL_G(N + OS2_G, :)         = tempL_G(N + OS2_G, :);

X_wH_G(Colorednodes_G{6}, :) = tempH_G(Colorednodes_G{6}, :);
X_wH_G(N + OS1_G, :)         = tempH_G(N + OS1_G, :);

% Synthesis
tempX_hatL_G = G0_G * X_wL_G;
tempX_hatH_G = G1_G * X_wH_G;

X_hat_G = tempX_hatL_G + tempX_hatH_G;

%% ============================================================
%% Step 5: Downsampling / redundancy removal
%% ============================================================

% Remove vertex-domain redundancy
X_hat_G_final = X_hat_G(1:N, :);
X_hat_G_final(connect_bpt_G, :) = ...
    0.5 * (X_hat_G(connect_bpt_G, :) + X_hat_G(N+1:end, :));

% Remove time-domain redundancy
X_hat_final = X_hat_G_final(:, 1:T);
X_hat_final(:, connect_bpt_T) = ...
    0.5 * (X_hat_G_final(:, connect_bpt_T) + X_hat_G_final(:, T+1:end));

% Final reconstructed signal
X_rec = X_hat_final;

%% ============================================================
%% Performance evaluation
%% ============================================================

f_hat = vec(X_rec);

MSE  = mean((f_hat_or - f_hat).^2);
NMSE = calculate_nmse(f_hat_or, f_hat);

P_rms = mean(f_hat_or.^2);
SNR   = 10 * log10(P_rms / MSE);

end
