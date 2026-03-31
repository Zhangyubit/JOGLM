function [X_rec, MSE, NMSE, SNR] = jtv_QMFfb_reconstruction(G, X)
% jtv_QMFfb_reconstruction - Joint time-vertex reconstruction using QMF graph filter banks
%
% This function performs joint time-vertex graph signal reconstruction using
% Quadrature Mirror Filter (QMF) banks defined on bipartite graphs.
%
% The reconstruction follows a separable structure:
%   1) Temporal-domain QMF filtering
%   2) Vertex-domain QMF filtering
%
% Inputs:
%   G     - Graph structure (GSPBox format)
%   X     - Graph signal matrix of size N x T
%
% Outputs:
%   X_rec - Reconstructed graph signal (N x T)
%   MSE   - Mean Squared Error
%   NMSE  - Normalized Mean Squared Error
%   SNR   - Signal-to-Noise Ratio (dB)
%
% Description:
%   - QMF filter banks use a single prototype kernel c(·), where:
%         Low-pass:  h(λ)
%         High-pass: h(2 - λ)
%   - This mirror structure ensures spectral complementarity on bipartite graphs.
%   - The algorithm applies time-domain processing followed by vertex-domain processing.

%% Parameter settings
filterlen = 12;

X_or = X;                      % Original signal
f_hat_or = vec(X_or);          % Vectorized reference signal

arange = [0 2];                % Spectral range

A_G = G.W;
N = G.N;
T = size(X,2);

% Construct temporal graph (ring graph)
GT = gsp_ring(T);
A_T = GT.W;

%% ============================================================
%% Step 1: Temporal domain QMF filtering
%% ============================================================

% Transpose for time-domain processing
X_T = X';   % Size: T x N

% Bipartite decomposition (time graph)
F_T = DSATUR(A_T);
[beta_T, bptG_T, beta_dist_T, Colorednodes_T] = harary_decomp(A_T, F_T);

theta_T = size(beta_T, 2);
Fmax_T  = size(beta_dist_T, 1);

% Construct QMF filter bank (single kernel)
[Ln_bpt_T, c_T] = QMFfb_filters(theta_T, bptG_T, filterlen, arange);

%% ---- Decomposition (time domain) ----
X_w_T = zeros(T, N, Fmax_T);
Channel_Name_T = cell(Fmax_T, 1);

for i = 1:Fmax_T
    if ~isempty(Colorednodes_T{i})
        tempX_w_T = X_T;
        
        for j = 1:theta_T
            if beta_dist_T(i,j) == 1
                % Low-pass filtering: h(λ)
                tempX_w_T = sgwt_cheby_op(tempX_w_T, Ln_bpt_T(:,:,j), c_T, arange);
                Channel_Name_T{i} = strcat(Channel_Name_T{i},'L');
            else
                % High-pass filtering: h(2 - λ)
                tempX_w_T = sgwt_cheby_op(tempX_w_T, 2*eye(T) - Ln_bpt_T(:,:,j), c_T, arange);
                Channel_Name_T{i} = strcat(Channel_Name_T{i},'H');
            end
        end
        
        % Sparse assignment
        X_w_T(Colorednodes_T{i}, :, i) = tempX_w_T(Colorednodes_T{i}, :, :);
    end
end

%% ---- Reconstruction (time domain) ----
X_hat_T = zeros(T, N, Fmax_T);

for i = 1:Fmax_T
    tempX_hat_T = X_w_T(:, :, i);
    
    for j = theta_T:-1:1
        if beta_dist_T(i,j) == 1
            tempX_hat_T = sgwt_cheby_op(tempX_hat_T, Ln_bpt_T(:,:,j), c_T, arange);
        else
            tempX_hat_T = sgwt_cheby_op(tempX_hat_T, 2*eye(T) - Ln_bpt_T(:,:,j), c_T, arange);
        end
    end
    
    X_hat_T(:,:,i) = tempX_hat_T;
end

% Aggregate channels and transpose back
X_rec_T = sum(X_hat_T,3)';   % Back to N x T


%% ============================================================
%% Step 2: Vertex domain QMF filtering
%% ============================================================

% Bipartite decomposition (vertex graph)
F_G = DSATUR(A_G);
[beta_G, bptG_G, beta_dist_G, Colorednodes_G] = harary_decomp(A_G, F_G);

theta_G = size(beta_G,2);
Fmax_G  = size(beta_dist_G,1);

% Construct QMF filter bank
[Ln_bpt_G, c_G] = QMFfb_filters(theta_G, bptG_G, filterlen, arange);

%% ---- Decomposition (vertex domain) ----
X_w_G = zeros(N, T, Fmax_G);
Channel_Name_G = cell(Fmax_G,1);

for i = 1:Fmax_G
    if ~isempty(Colorednodes_G{i})
        tempX_w_G = X_rec_T;
        
        for j = 1:theta_G
            if beta_dist_G(i,j) == 1
                % Low-pass
                tempX_w_G = sgwt_cheby_op(tempX_w_G, Ln_bpt_G(:,:,j), c_G, arange);
                Channel_Name_G{i} = strcat(Channel_Name_G{i},'L');
            else
                % High-pass (mirror)
                tempX_w_G = sgwt_cheby_op(tempX_w_G, 2*eye(N) - Ln_bpt_G(:,:,j), c_G, arange);
                Channel_Name_G{i} = strcat(Channel_Name_G{i},'H');
            end
        end
        
        % Sparse assignment
        X_w_G(Colorednodes_G{i}, :, i) = tempX_w_G(Colorednodes_G{i}, :, :);
    end
end

%% ---- Reconstruction (vertex domain) ----
X_hat_G = zeros(N, T, Fmax_G);

for i = 1:Fmax_G
    tempX_hat_G = X_w_G(:, :, i);
    
    for j = theta_G:-1:1
        if beta_dist_G(i,j) == 1
            tempX_hat_G = sgwt_cheby_op(tempX_hat_G, Ln_bpt_G(:,:,j), c_G, arange);
        else
            tempX_hat_G = sgwt_cheby_op(tempX_hat_G, 2*eye(N) - Ln_bpt_G(:,:,j), c_G, arange);
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