%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Copyright (c) 2025 Yu Zhang
%
%  Description:
%  This file denoises an image sequence using the joint time-vertex 
%  oversampled graph Laplacian framework.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear;

gsp_start

%% Parameter settings
sigma = 10;                % Noise standard deviation
filterlen = 12;            % Filter length
norm_type = 'sym';         % Type of Laplacian normalization

%% Step 1: Construct graph structure
% Use a 2D grid graph to model the spatial structure of the image
N = 35;                    % Image size (N x N grid)
G = gsp_2dgrid(N);         % Generate 2D grid graph

% Obtain graph Laplacian
L = G.L;                   % Graph Laplacian matrix
% Alternatively: A = G.W;

%% Step 2: Graph Fourier Transform (GFT)
% Compute eigenvectors and eigenvalues of the graph Laplacian
[GFT_G, V_G, Lambda_G] = gft(L);

%% Step 3: Time domain (FRFT-like construction)
T = 101;                   % Number of time samples (interpolation points)
x = linspace(-1, 1, T);    % Parameter domain (acts like time/FRFT domain)

GT = gsp_ring(T);          % Construct ring graph for time domain
[GFT_T, V_T, Lambda_T] = gft(GT.L);

%% Step 4: Define signals (image interpolation)
% Load two images (e.g., digit "0" and digit "6")
f0 = imread('digit0.png'); 
f6 = imread('digit6.png');

% Resize and convert to double precision
f0 = double(imresize(f0, [N N]));
f6 = double(imresize(f6, [N N]));

% Store original clean signals
f0_original = f0;
f6_original = f6;

%% Construct original (clean) joint signal
% Linear interpolation between f0 and f6
f_original = @(x) f0_original + (f6_original - f0_original) * (x + 1) / 2;

% Generate signal matrix: each column corresponds to one "time" sample
F_x_original = arrayfun(@(t) reshape(f_original(t), [], 1), x, 'UniformOutput', false);
signal_matrix_original = cell2mat(F_x_original); % Size: (N^2 × T)

% Vectorized clean signal
f_or = vec(signal_matrix_original);

%% Add Gaussian noise
rng(42);                              % Fix random seed for reproducibility
f0 = f0 + sigma * randn(size(f0));    % Add noise to first image
f6 = f6 + sigma * randn(size(f6));    % Add noise to second image

% Define noisy interpolation function
f = @(x) f0 + (f6 - f0) * (x + 1) / 2;

% Generate noisy signal matrix
F_x = arrayfun(@(t) reshape(f(t), [], 1), x, 'UniformOutput', false);
signal_matrix = cell2mat(F_x);

% Vectorized noisy signal
f_noise = vec(signal_matrix);

%% Step 5: Denoising using graph filter banks

% Oversampled graph filter bank denoising
X_rec_os = jtv_ogfb_denoise(G, signal_matrix, filterlen, norm_type, sigma);

% Critically sampled graph filter bank denoising
X_rec_or = jtv_gfb_denoise(G, signal_matrix, filterlen, norm_type, sigma);

%% Sanity check on output dimensions
n = N * N;
T = length(x);

if size(X_rec_os,1) ~= n || size(X_rec_os,2) ~= T
    error('X_rec_os has incorrect size, expected [N^2, T]');
end
if size(X_rec_or,1) ~= n || size(X_rec_or,2) ~= T
    error('X_rec_or has incorrect size, expected [N^2, T]');
end

%% Step 6: Vectorization (column-major order)
f_hat_os = X_rec_os(:);    % Oversampled reconstruction
f_hat_or = X_rec_or(:);    % Critically sampled reconstruction

%% Step 7: Compute MSE and SNR
MSE_noisy = mean((f_or - f_noise).^2);
SNR_noisy = 10 * log10(mean(f_or.^2) / MSE_noisy);

MSE_os = mean((f_or - f_hat_os).^2);
SNR_os = 10 * log10(mean(f_or.^2) / MSE_os);

MSE_or = mean((f_or - f_hat_or).^2);
SNR_or = 10 * log10(mean(f_or.^2) / MSE_or);

% Display results
disp('1 to 7')
fprintf('sigma value = %s\n', rats(sigma));
fprintf('MSE of noisy signal = %.4f, SNR = %.2f dB\n', MSE_noisy, SNR_noisy);
fprintf('MSE of denoised signal by Oversampling = %.4f, SNR = %.2f dB\n', MSE_os, SNR_os);
fprintf('MSE of denoised signal by Critically Sampling = %.4f, SNR = %.2f dB\n', MSE_or, SNR_or);

%% Step 8: Recover specific frames (x = -1 and x = +1)
[~, idx_minus1] = min(abs(x - (-1)));
[~, idx_plus1 ] = min(abs(x -   1));

% Oversampled reconstruction
img_os_minus1 = reshape(X_rec_os(:, idx_minus1), [N, N]);
img_os_plus1  = reshape(X_rec_os(:, idx_plus1),  [N, N]);

% Critically sampled reconstruction
img_or_minus1 = reshape(X_rec_or(:, idx_minus1), [N, N]);
img_or_plus1  = reshape(X_rec_or(:, idx_plus1),  [N, N]);

%% Step 9: Visualization
figure;
subplot(2,4,1); imshow(f0_original,[]); title('orig f0');
subplot(2,4,2); imshow(f0,[]);          title('noisy f0');
subplot(2,4,3); imshow(img_or_minus1,[]); title('OR recon f0');
subplot(2,4,4); imshow(img_os_minus1,[]); title('OS recon f0');

subplot(2,4,5); imshow(f6_original,[]); title('orig f6');
subplot(2,4,6); imshow(f6,[]);          title('noisy f6');
subplot(2,4,7); imshow(img_or_plus1,[]); title('OR recon f6');
subplot(2,4,8); imshow(img_os_plus1,[]); title('OS recon f6');

%% Step 10: Save results as images

% Case: digits 0–6
imwrite(mat2gray(f0_original), 'orig_f0.jpg');
imwrite(mat2gray(f0),          'noisy_f0.jpg');
imwrite(mat2gray(img_or_minus1), 'OR_recon_f0.jpg');
imwrite(mat2gray(img_os_minus1), 'OS_recon_f0.jpg');

imwrite(mat2gray(f6_original), 'orig_f6.jpg');
imwrite(mat2gray(f6),          'noisy_f6.jpg');
imwrite(mat2gray(img_or_plus1), 'OR_recon_f6.jpg');
imwrite(mat2gray(img_os_plus1), 'OS_recon_f6.jpg');

% Case: digits 1–7 (renamed outputs)
imwrite(mat2gray(f0_original), 'orig_f1.jpg');
imwrite(mat2gray(f0),          'noisy_f1.jpg');
imwrite(mat2gray(img_or_minus1), 'OR_recon_f1.jpg');
imwrite(mat2gray(img_os_minus1), 'OS_recon_f1.jpg');

imwrite(mat2gray(f6_original), 'orig_f7.jpg');
imwrite(mat2gray(f6),          'noisy_f7.jpg');
imwrite(mat2gray(img_or_plus1), 'OR_recon_f7.jpg');
imwrite(mat2gray(img_os_plus1), 'OS_recon_f7.jpg');
