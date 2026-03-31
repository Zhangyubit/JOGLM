function [X_rec, MSE, NMSE, SNR] = jtv_Bifb_reconstruction(G, X)
% jtv_Bifb_reconstruction - Joint time-vertex reconstruction using biorthogonal graph filter banks
%
% This function performs joint time-vertex graph signal reconstruction by
% applying a separable two-stage processing:
%   1) Temporal graph filter bank (on time domain)
%   2) Vertex-domain graph filter bank
%
% Both stages use biorthogonal graph filter banks constructed on bipartite graphs.
%
% Input:
%   G     - Graph structure (GSPBox format)
%   X     - Input graph signal matrix of size N x T
%
% Outputs:
%   X_rec - Reconstructed graph signal (N x T)
%   MSE   - Mean Squared Error
%   NMSE  - Normalized Mean Squared Error
%   SNR   - Signal-to-Noise Ratio (dB)
%
% Description:
%   - The algorithm first processes the signal along the time dimension
%     (treating each node as a temporal signal), and then along the vertex
%     dimension (treating each time instant as a graph signal).
%   - Each stage consists of:
%         Decomposition → Channel-wise processing → Reconstruction
%   - Perfect reconstruction is enabled via biorthogonal filter design.

%% Parameter settings
filterlen = 12;
X_or = X;                          % Original signal
f_hat_or = vec(X_or);              % Vectorized reference signal

arange = [0 2];                    % Spectral range

A_G = G.W;
N = G.N;
T = size(X,2);

% Construct temporal graph (ring graph)
GT = gsp_ring(T);
A_T = GT.W;

%% ============================================================
%% Step 1: Temporal domain processing (decomposition → reconstruction)
%% ============================================================

% Transpose signal: process along time for each node
X_T = X';   % Now size is T x N

% Graph coloring and bipartite decomposition (time graph)
F_T = DSATUR(A_T);
[beta_T, bptG_T, beta_dist_T, Colorednodes_T] = harary_decomp(A_T, F_T);

theta_T = size(beta_T, 2);         % Number of decomposition stages
Fmax_T  = size(beta_dist_T, 1);    % Number of channels

% Construct biorthogonal filter bank
[Ln_bpt_T, c_d_T, c_r_T] = Bifb_filters(theta_T, bptG_T, filterlen, arange);

%% ---- Decomposition (time domain) ----
X_w_T = zeros(T, N, Fmax_T);       % Channel outputs
Channel_Name_T = cell(Fmax_T, 1);

for i = 1:Fmax_T
    if ~isempty(Colorednodes_T{i})
        tempX_w_T = X_T;
        
        for j = 1:theta_T
            if beta_dist_T(i, j) == 1
                % Low-pass filtering
                tempX_w_T = sgwt_cheby_op(tempX_w_T, Ln_bpt_T(:,:,j), c_d_T{1}, arange);
                Channel_Name_T{i} = strcat(Channel_Name_T{i}, 'L');
            else
                % High-pass filtering
                tempX_w_T = sgwt_cheby_op(tempX_w_T, Ln_bpt_T(:,:,j), c_d_T{2}, arange);
                Channel_Name_T{i} = strcat(Channel_Name_T{i}, 'H');
            end
        end
        
        % Store only selected nodes (sparse assignment)
        X_w_T(Colorednodes_T{i}, :, i) = tempX_w_T(Colorednodes_T{i}, :, :);
    end
end

%% ---- Reconstruction (time domain) ----
X_hat_T = zeros(T, N, Fmax_T);

for i = 1:Fmax_T
    tempX_hat_T = X_w_T(:, :, i);
    
    for j = theta_T:-1:1
        if beta_dist_T(i, j) == 1
            tempX_hat_T = sgwt_cheby_op(tempX_hat_T, Ln_bpt_T(:,:,j), c_r_T{1}, arange);
        else
            tempX_hat_T = sgwt_cheby_op(tempX_hat_T, Ln_bpt_T(:,:,j), c_r_T{2}, arange);
        end
    end
    
    X_hat_T(:,:,i) = tempX_hat_T;
end

% Aggregate all channels and transpose back
X_rec_T = sum(X_hat_T,3)';   % Back to N x T


%% ============================================================
%% Step 2: Vertex domain processing (decomposition → reconstruction)
%% ============================================================

% Graph coloring and bipartite decomposition (vertex graph)
F_G = DSATUR(A_G);
[beta_G, bptG_G, beta_dist_G, Colorednodes_G] = harary_decomp(A_G, F_G);

theta_G = size(beta_G,2);
Fmax_G  = size(beta_dist_G,1);

% Construct filter bank
[Ln_bpt_G, c_d_G, c_r_G] = Bifb_filters(theta_G, bptG_G, filterlen, arange);

%% ---- Decomposition (vertex domain) ----
X_w_G = zeros(N, T, Fmax_G);
Channel_Name_G = cell(Fmax_G,1);

for i = 1:Fmax_G
    if ~isempty(Colorednodes_G{i})
        tempX_w_G = X_rec_T;
        
        for j = 1:theta_G
            if beta_dist_G(i,j) == 1
                tempX_w_G = sgwt_cheby_op(tempX_w_G, Ln_bpt_G(:,:,j), c_d_G{1}, arange);
                Channel_Name_G{i} = strcat(Channel_Name_G{i},'L');
            else
                tempX_w_G = sgwt_cheby_op(tempX_w_G, Ln_bpt_G(:,:,j), c_d_G{2}, arange);
                Channel_Name_G{i} = strcat(Channel_Name_G{i},'H');
            end
        end
        
        % Sparse assignment
        X_w_G(Colorednodes_G{i}, :, i) = tempX_w_G(Colorednodes_G{i},:,:);
    end
end

%% ---- Reconstruction (vertex domain) ----
X_hat_G = zeros(N, T, Fmax_G);

for i = 1:Fmax_G
    tempX_hat_G = X_w_G(:, :, i);
    
    for j = theta_G:-1:1
        if beta_dist_G(i,j) == 1
           tempX_hat_G = sgwt_cheby_op(tempX_hat_G, Ln_bpt_G(:,:,j), c_r_G{1}, arange);
        else
           tempX_hat_G = sgwt_cheby_op(tempX_hat_G, Ln_bpt_G(:,:,j), c_r_G{2}, arange);
        end
    end
    
    X_hat_G(:,:,i) = tempX_hat_G;
end

%% ============================================================
%% Final reconstruction
%% ============================================================

X_rec = sum(X_hat_G,3);

%% ============================================================
%% Performance evaluation
%% ============================================================

f_hat = vec(X_rec);

MSE  = mean((f_hat_or - f_hat).^2);
NMSE = calculate_nmse(f_hat_or, f_hat);

P_rms = mean(f_hat_or.^2);
SNR   = 10 * log10(P_rms / MSE);

end