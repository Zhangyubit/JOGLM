clear; clc; close all;

disp('***** Denoising with Oversampled Graph Filter Bank and other methods *****')

%% ============================================================
%% Parameter settings
%% ============================================================

sigma     = 1;        % Noise standard deviation
filterlen = 12;       % Filter length
norm_type = 'sym';    % Laplacian normalization type
thresh    = 3;        % Threshold for denoising

%% ============================================================
%% Section 1: Graph construction and signal generation
%% ============================================================

gsp_start;

% Load graph structure
% G = gsp_minnesota();   % Example graph
G = gsp_Ex();            % Custom example graph

A  = G.W;
xy = G.coords;
N  = G.N;

% Load original graph signal
% load min_graph_signal
load Ex_signal6   % Provides signal f

T = 3;   % Number of time samples. Notes: or T = 4; 

% Construct time-extended signal (same signal replicated across time)
X = repmat(f, 1, T);

% Store original signal for evaluation
X_or = X;
f_or = vec(X_or);

%% Add Gaussian noise

rng(42);   % Fix random seed for reproducibility
X = X + sigma * randn(size(X));

f_noise = vec(X);

%% ============================================================
%% Section 2: Graph filter bank denoising
%% ============================================================

% Apply different graph filter bank methods

% Oversampled Graph Filter Bank (OGFB)
X_rec_os = jtv_ogfb_denoise(G, X, filterlen, norm_type, sigma, thresh);

% Critically-sampled Graph Filter Bank (GFB)
X_rec_or = jtv_gfb_denoise(G, X, filterlen, norm_type, sigma, thresh);

% Biorthogonal Graph Filter Bank (BiFB)
X_rec_bi = jtv_Bifb_denoise(G, X, filterlen, sigma, thresh);

% Quadrature Mirror Filter Bank (QMF)
X_rec_QMF = jtv_QMFfb_denoise(G, X, filterlen, sigma, thresh);

%% ============================================================
%% Section 3: Performance evaluation
%% ============================================================

% Vectorize reconstructed signals
f_hat_os  = vec(X_rec_os);
f_hat_or  = vec(X_rec_or);
f_hat_bi  = vec(X_rec_bi);
f_hat_QMF = vec(X_rec_QMF);

% Compute MSE and SNR
P_rms = mean(f_or.^2);

MSE_os  = mean((f_or - f_hat_os).^2);
SNR_os  = 10 * log10(P_rms / MSE_os);

MSE_or  = mean((f_or - f_hat_or).^2);
SNR_or  = 10 * log10(P_rms / MSE_or);

MSE_bi  = mean((f_or - f_hat_bi).^2);
SNR_bi  = 10 * log10(P_rms / MSE_bi);

MSE_QMF = mean((f_or - f_hat_QMF).^2);
SNR_QMF = 10 * log10(P_rms / MSE_QMF);

% Noisy signal baseline
MSE_noisy = mean((f_or - f_noise).^2);
SNR_noisy = 10 * log10(P_rms / MSE_noisy);

%% ============================================================
%% Display results
%% ============================================================

disp('***** Denoising Performance *****')

fprintf('Noise level (sigma) = %s\n', rats(sigma));

fprintf('MSE of noisy signal = %.4f\n', MSE_noisy);
fprintf('SNR of noisy signal = %.2f dB\n\n', SNR_noisy);

fprintf('MSE (Biorthogonal GFB) = %.4f\n', MSE_bi);
fprintf('SNR (Biorthogonal GFB) = %.2f dB\n\n', SNR_bi);

fprintf('MSE (QMF GFB) = %.4f\n', MSE_QMF);
fprintf('SNR (QMF GFB) = %.2f dB\n\n', SNR_QMF);

fprintf('MSE (Critically-sampled GFB) = %.4f\n', MSE_or);
fprintf('SNR (Critically-sampled GFB) = %.2f dB\n\n', SNR_or);

fprintf('MSE (Oversampled GFB) = %.4f\n', MSE_os);
fprintf('SNR (Oversampled GFB) = %.2f dB\n', SNR_os);

%% ============================================================
%% Section 4: Visualization (optional)
%% ============================================================
% Note: Visualization may be time-consuming

plot_prod_graph_signal(G, T, f_or);        % Original signal
plot_prod_graph_signal(G, T, f_noise);     % Noisy signal
plot_prod_graph_signal(G, T, f_hat_QMF);   % QMF result
plot_prod_graph_signal(G, T, f_hat_bi);    % Biorthogonal result
plot_prod_graph_signal(G, T, f_hat_or);    % Critically sampled result
plot_prod_graph_signal(G, T, f_hat_os);    % Oversampled result