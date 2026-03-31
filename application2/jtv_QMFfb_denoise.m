function X_rec = jtv_QMFfb_denoise(G, X, filterlen, sigma, thresh)
disp('***** Joint Time-Vertex QMF Filter Bank Denoising *****')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameter settings
arange = [0 2];                 % Spectral approximation range

A_G = G.W;                      % Adjacency matrix of the vertex graph
N = G.N;                        % Number of vertices

T = size(X,2);                  % Number of time samples

GT = gsp_ring(T);               % Construct time graph (ring structure)
A_T = GT.W;                     % Adjacency matrix of time graph

% Nonlinear scaling of threshold (empirical strategy)
ord = 0.7;
thresh = thresh.^ord;
thresh = min(thresh, 3^ord);    % Upper bound for stability

%% ================= Filterbank Implementation =================

%% ---- Time-domain decomposition - thresholding - reconstruction ----
X_T = X';  % Transpose: now T × N, each column is a temporal signal on a vertex

% Construct bipartite graph filter banks on the time graph
F_T = DSATUR(A_T);   % Graph coloring
[beta_T, bptG_T, beta_dist_T, Colorednodes_T] = harary_decomp(A_T, F_T);

theta_T = size(beta_T, 2);      % Number of decomposition stages
Fmax_T = size(beta_dist_T, 1);  % Number of output channels

[Ln_bpt_T, c_T] = QMFfb_filters(theta_T, bptG_T, filterlen, arange);

% Initialize filterbank outputs
X_w_T = zeros(T, N, Fmax_T);
Channel_Name_T = cell(Fmax_T, 1);

% Decomposition: apply filters along time dimension
for i = 1:Fmax_T
    if ~isempty(Colorednodes_T{i})
        tempX_w_T = X_T;
        for j = 1:theta_T
            if beta_dist_T(i,j) == 1
                % Low-pass branch
                tempX_w_T = sgwt_cheby_op(tempX_w_T, Ln_bpt_T(:,:,j), c_T, arange);
                Channel_Name_T{i} = strcat(Channel_Name_T{i}, 'L');
            else
                % High-pass branch (QMF structure)
                tempX_w_T = sgwt_cheby_op(tempX_w_T, 2*eye(T) - Ln_bpt_T(:,:,j), c_T, arange);
                Channel_Name_T{i} = strcat(Channel_Name_T{i}, 'H');
            end
        end
        % Sparse assignment based on coloring
        X_w_T(Colorednodes_T{i},:,i) = tempX_w_T(Colorednodes_T{i},:,:);
    end
end

% Thresholding: applied only to high-frequency channels
thresh_T = thresh * sigma; 
for i = 1:Fmax_T
    if contains(char(Channel_Name_T{i}), 'H')
        X_w_T(:,:,i) = X_w_T(:,:,i) .* (abs(X_w_T(:,:,i)) > thresh_T);
    end
end

%% ---- Time-domain reconstruction ----
X_hat_T = zeros(T, N, Fmax_T);

for i = 1:Fmax_T
    tempX_hat_T = X_w_T(:,:,i);
    for j = theta_T:-1:1
        if beta_dist_T(i,j) == 1
            tempX_hat_T = sgwt_cheby_op(tempX_hat_T, Ln_bpt_T(:,:,j), c_T, arange);
        else
            tempX_hat_T = sgwt_cheby_op(tempX_hat_T, 2*eye(T) - Ln_bpt_T(:,:,j), c_T, arange);
        end
    end
    X_hat_T(:,:,i) = tempX_hat_T;
end

% Aggregate channels and transpose back to N × T
X_rec_T = sum(X_hat_T,3)';

%% ---- Vertex-domain decomposition - thresholding - reconstruction ----
F_G = DSATUR(A_G);
[beta_G, bptG_G, beta_dist_G, Colorednodes_G] = harary_decomp(A_G, F_G);

theta_G = size(beta_G,2);
Fmax_G = size(beta_dist_G,1);

[Ln_bpt_G, c_G] = QMFfb_filters(theta_G, bptG_G, filterlen, arange);

% Decomposition
X_w_G = zeros(N, T, Fmax_G);
Channel_Name_G = cell(Fmax_G,1);

for i = 1:Fmax_G
    if ~isempty(Colorednodes_G{i})
        tempX_w_G = X_rec_T;
        for j = 1:theta_G
            if beta_dist_G(i,j) == 1
                tempX_w_G = sgwt_cheby_op(tempX_w_G, Ln_bpt_G(:,:,j), c_G, arange);
                Channel_Name_G{i} = strcat(Channel_Name_G{i}, 'L');
            else
                tempX_w_G = sgwt_cheby_op(tempX_w_G, 2*eye(N) - Ln_bpt_G(:,:,j), c_G, arange);
                Channel_Name_G{i} = strcat(Channel_Name_G{i}, 'H');
            end
        end
        X_w_G(Colorednodes_G{i},:,i) = tempX_w_G(Colorednodes_G{i},:,:);
    end
end

% Thresholding (high-frequency components only)
thresh_G = thresh * sigma; 
for i = 1:Fmax_G
    if contains(char(Channel_Name_G{i}), 'H')
        X_w_G(:,:,i) = X_w_G(:,:,i) .* (abs(X_w_G(:,:,i)) > thresh_G);
    end
end

%% ---- Vertex-domain reconstruction ----
X_hat_G = zeros(N, T, Fmax_G);

for i = 1:Fmax_G
    tempX_hat_G = X_w_G(:,:,i);
    for j = theta_G:-1:1
        if beta_dist_G(i,j) == 1
            tempX_hat_G = sgwt_cheby_op(tempX_hat_G, Ln_bpt_G(:,:,j), c_G, arange);
        else
            tempX_hat_G = sgwt_cheby_op(tempX_hat_G, 2*eye(N) - Ln_bpt_G(:,:,j), c_G, arange);
        end
    end
    X_hat_G(:,:,i) = tempX_hat_G;
end

%% ================= Final Reconstruction =================
% Sum over all channels to obtain final denoised signal
X_rec = sum(X_hat_G,3);
end