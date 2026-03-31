function X_rec = jtv_Bifb_denoise(G, X, filterlen, sigma, thresh)
% jtv_Bifb_denoise - Joint time-vertex denoising using biorthogonal graph filter banks
%
% This function performs denoising of a time-vertex graph signal using a 
% two-stage (time → vertex) biorthogonal graph filter bank (BiFB) framework.
% Thresholding is applied to high-frequency components in both domains.
%
% Input:
%   G          - Graph structure (from GSPBox)
%   X          - Graph signal matrix of size N × T
%   filterlen  - Length of the graph filters
%   sigma      - Noise standard deviation
%   thresh     - Threshold scaling factor
%
% Output:
%   X_rec      - Reconstructed (denoised) graph signal

disp('***** Joint Time-Vertex Biorthogonal Filter Bank Denoising *****')

%% Parameters
arange = [0 2];            % Spectral approximation range

A_G = G.W;                 % Adjacency matrix (vertex domain)
N = G.N;                   % Number of vertices
T = size(X,2);             % Number of time samples

% Construct time graph (ring structure)
GT = gsp_ring(T);
A_T = GT.W;


%% ================= TIME DOMAIN PROCESSING =================
% Decomposition → Thresholding → Reconstruction

% Transpose signal: now each column is a temporal signal on a vertex
X_T = X';                  % Size: T × N

% Graph coloring and Harary decomposition
F_T = DSATUR(A_T);
[beta_T, bptG_T, beta_dist_T, Colorednodes_T] = harary_decomp(A_T, F_T);

theta_T = size(beta_T, 2);         % Number of decomposition stages
Fmax_T  = size(beta_dist_T, 1);    % Number of channels

% Construct biorthogonal filter bank
[Ln_bpt_T, c_d_T, c_r_T] = Bifb_filters(theta_T, bptG_T, filterlen, arange);

% Initialize filter bank outputs
X_w_T = zeros(T, N, Fmax_T);
Channel_Name_T = cell(Fmax_T, 1);

%% ----- Decomposition (along time for each vertex) -----
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
        % Sparse assignment (only write to selected nodes)
        X_w_T(Colorednodes_T{i}, :, i) = tempX_w_T(Colorednodes_T{i}, :, :);
    end
end

%% ----- Thresholding (time domain high-frequency components) -----
thresh_T = thresh * sigma;

for i = 1:Fmax_T
    if ~isempty(Channel_Name_T{i}) && contains(char(Channel_Name_T{i}), 'H')
        X_w_T(:, :, i) = X_w_T(:, :, i) .* ...
                         (abs(X_w_T(:, :, i)) > thresh_T);
    end
end

%% ----- Reconstruction (time domain) -----
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

% Aggregate channels and transpose back
X_rec_T = sum(X_hat_T,3)';   % T × N → N × T


%% ================= VERTEX DOMAIN PROCESSING =================
% Decomposition → Thresholding → Reconstruction

F_G = DSATUR(A_G);
[beta_G, bptG_G, beta_dist_G, Colorednodes_G] = harary_decomp(A_G, F_G);

theta_G = size(beta_G,2);
Fmax_G  = size(beta_dist_G,1);

[Ln_bpt_G, c_d_G, c_r_G] = Bifb_filters(theta_G, bptG_G, filterlen, arange);

% Initialize outputs
X_w_G = zeros(N, T, Fmax_G);
Channel_Name_G = cell(Fmax_G,1);

%% ----- Decomposition (vertex domain) -----
for i = 1:Fmax_G
    if ~isempty(Colorednodes_G{i})
        tempX_w_G = X_rec_T;
        for j = 1:theta_G
            if beta_dist_G(i,j) == 1
                tempX_w_G = sgwt_cheby_op(tempX_w_G, Ln_bpt_G(:,:,j), c_d_G{1}, arange);
                Channel_Name_G{i} = strcat(Channel_Name_G{i}, 'L');
            else
                tempX_w_G = sgwt_cheby_op(tempX_w_G, Ln_bpt_G(:,:,j), c_d_G{2}, arange);
                Channel_Name_G{i} = strcat(Channel_Name_G{i}, 'H');
            end
        end
        X_w_G(Colorednodes_G{i}, :, i) = tempX_w_G(Colorednodes_G{i}, :, :);
    end
end

%% ----- Thresholding (vertex domain high-frequency components) -----
thresh_G = thresh * sigma;

for i = 1:Fmax_G
    if ~isempty(Channel_Name_G{i}) && contains(char(Channel_Name_G{i}), 'H')
        X_w_G(:, :, i) = X_w_G(:, :, i) .* ...
                         (abs(X_w_G(:, :, i)) > thresh_G);
    end
end

%% ----- Reconstruction (vertex domain) -----
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

%% ================= FINAL RECONSTRUCTION =================
% Sum over all channels
X_rec = sum(X_hat_G,3);

end