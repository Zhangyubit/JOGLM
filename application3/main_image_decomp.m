%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Copyright (c) 2025 Yu Zhang
%
%  Description:
%  This script performs image decomposition using a joint time-vertex
%  oversampled graph Laplacian framework. The input images are interpolated
%  to form a joint signal, which is then processed via graph filter banks.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear;

gsp_start

%% ================= Step 1: Construct Graph Structure =================
% Use a 2D grid graph to represent image pixels
N = 28;                         % Image size: 28 × 28
G = gsp_2dgrid(N);              % Generate 2D grid graph

% Graph Laplacian
L = G.L;                        % Graph Laplacian matrix

%% ================= Step 2: Graph Fourier Transform =================
% Compute eigenvectors and eigenvalues of the Laplacian
[GFT_G, V_G, Lambda_G] = gft(L);

%% ================= Step 3: Time (FRFT-like) Domain Setup =================
T = 100;                        % Number of interpolation points
x = linspace(-1, 1, T);         % Parameter domain

GT = gsp_ring(T);               % Construct time graph (ring structure)
[GFT_T, V_T, Lambda_T] = gft(GT.L);

%% ================= Step 4: Define Signals =================
% Load two images (digit "0" and "6")
f0 = imread('digit0.png');
f6 = imread('digit6.png');

% Resize and convert to double precision
f0 = double(imresize(f0, [N N]));
f6 = double(imresize(f6, [N N]));

% Store original clean images
f0_original = f0;
f6_original = f6;

%% ================= Step 5: Add Noise =================
sigma = 5;     
rng(42);                               % Fix random seed

f0 = f0 + sigma * randn(size(f0));      % Add Gaussian noise
f6 = f6 + sigma * randn(size(f6));

%% ================= Step 6: Construct Interpolated Signal =================
% Define interpolation function between f0 and f6
f = @(x) f0 + (f6 - f0) * (x + 1) / 2;

% Generate signals for each interpolation parameter x
F_x = arrayfun(@(t) reshape(f(t), [], 1), x, 'UniformOutput', false);

% Form joint signal matrix (each column corresponds to one x)
signal_matrix = cell2mat(F_x);

%% ================= Step 7: Reconstruction via Oversampled GFB =================
[hatX_os, X_rec_os, MSE_os, NMSE_os, SNR_os] = ...
    jtv_ogfb_reconstruction(G, signal_matrix);

%% ================= Step 8: Dimension Check =================
n = N * N;
[nt_rows, nt_cols] = size(X_rec_os);

if nt_rows ~= n || nt_cols ~= length(x)
    error('X_rec_os has incorrect size: expected [N^2, T]');
end

%% ================= Step 9: Vectorization =================
% Column-major vectorization (MATLAB default)
f_colmajor = X_rec_os(:);   
% Interpretation: pixels are stacked frame by frame (time varies slowest)

% Time-fast ordering (time index varies fastest)
f_timefast = reshape(X_rec_os.', [], 1);

%% ================= Step 10: Recover Individual Frames =================
j = 1; % Example: first frame
img_j = reshape(X_rec_os(:, j), [N, N]);
img_j = real(img_j);  % Remove numerical imaginary parts if any

%% ================= Step 11: Recover Boundary Frames =================
% Find indices corresponding to x = -1 and x = +1
[~, idx_minus1] = min(abs(x - (-1)));
[~, idx_plus1 ] = min(abs(x -   1));

img_minus1 = reshape(X_rec_os(:, idx_minus1), [N, N]);
img_plus1  = reshape(X_rec_os(:, idx_plus1),  [N, N]);

img_minus1 = real(img_minus1);
img_plus1  = real(img_plus1);

% Compute reconstruction errors (sanity check)
err0 = max(abs(img_minus1(:) - f0_original(:)));
err6 = max(abs(img_plus1(:)  - f6_original(:)));

fprintf('max error for x=-1: %g, for x=+1: %g\n', err0, err6);

%% ================= Step 12: Construct 3D Image Volume =================
imgVol = zeros(N, N, length(x));

for t = 1:length(x)
    imgVol(:, :, t) = reshape(X_rec_os(:, t), [N, N]);
end

% Equivalent one-line version:
% imgVol = reshape(X_rec_os, [N, N, length(x)]);

%% ================= Step 13: Visualization =================
figure;

subplot(2,3,1); imshow(f0_original,[]); title('Original f0');
subplot(2,3,2); imshow(f0,[]);          title('Noisy f0');
subplot(2,3,3); imshow(img_minus1,[]);  
title(sprintf('Reconstruction at x = -1 (idx = %d)', idx_minus1));

subplot(2,3,4); imshow(f6_original,[]); title('Original f6');
subplot(2,3,5); imshow(f6,[]);          title('Noisy f6');
subplot(2,3,6); imshow(img_plus1,[]);   
title(sprintf('Reconstruction at x = +1 (idx = %d)', idx_plus1));