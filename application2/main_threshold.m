clear; clc; close all;

disp('***** Denoising with different thresholds *****')

%% ============================================================
%% Parameter settings
%% ============================================================

sigma       = 1;          % Noise standard deviation
filterlen   = 12;         % Filter length
norm_type   = 'sym';      % Laplacian normalization type
thresh_values = 0:0.2:4;  % Range of threshold values

%% ============================================================
%% Section 1: Graph construction
%% ============================================================

gsp_start;

% Load graph structure
G = gsp_Ex();
% G = gsp_minnesota();   % Alternative graph

A  = G.W;
xy = G.coords;
N  = G.N;

%% ============================================================
%% Section 2: Signal construction
%% ============================================================

% Load graph signal
load Ex_signal6
% load min_graph_signal

T = 3;   % Number of time samples

% Replicate the signal across time
X = repmat(f, 1, T);

% Store original signal
X_or = X;
f_or = vec(X_or);

%% ============================================================
%% Section 3: Add Gaussian noise
%% ============================================================

rng(42);   % Fix random seed for reproducibility
X = X + sigma * randn(size(X));

%% ============================================================
%% Section 4: Threshold sweep (performance evaluation)
%% ============================================================

num_thresh = length(thresh_values);

% Pre-allocate MSE storage
MSE_os  = zeros(num_thresh,1);
MSE_or  = zeros(num_thresh,1);
MSE_bi  = zeros(num_thresh,1);
MSE_QMF = zeros(num_thresh,1);

% Loop over threshold values
for i = 1:num_thresh
    
    thresh = thresh_values(i);
    
    % ----- Oversampled Graph Filter Bank -----
    X_rec_os = jtv_ogfb_denoise(G, X, filterlen, norm_type, sigma, thresh);
    f_hat_os = vec(X_rec_os);
    MSE_os(i) = mean((f_or - f_hat_os).^2);
    
    % ----- Critically sampled Graph Filter Bank -----
    X_rec_or = jtv_gfb_denoise(G, X, filterlen, norm_type, sigma, thresh);
    f_hat_or = vec(X_rec_or);
    MSE_or(i) = mean((f_or - f_hat_or).^2);
    
    % ----- Biorthogonal Graph Filter Bank -----
    X_rec_bi = jtv_Bifb_denoise(G, X, filterlen, sigma, thresh);
    f_hat_bi = vec(X_rec_bi);
    MSE_bi(i) = mean((f_or - f_hat_bi).^2);
    
    % ----- QMF-based Graph Filter Bank -----
    X_rec_QMF = jtv_QMFfb_denoise(G, X, filterlen, sigma, thresh);
    f_hat_QMF = vec(X_rec_QMF);
    MSE_QMF(i) = mean((f_or - f_hat_QMF).^2);
    
end

%% ============================================================
%% Section 5: Display results
%% ============================================================

disp('***** MSE for different thresholds *****')

TBL = table(thresh_values', MSE_or, MSE_os, MSE_bi, MSE_QMF, ...
    'VariableNames', {'Threshold','CriticallySampled',...
                      'Oversampled','GraphBior','GraphQMF'});

disp(TBL)

%% ============================================================
%% Section 6: Visualization
%% ============================================================

figure('Position', [200, 200, 650, 650]); % Enlarged figure size

plot(thresh_values, MSE_QMF,  '-o', 'LineWidth', 4, 'MarkerSize', 12); hold on;
plot(thresh_values, MSE_bi,   '-s', 'LineWidth', 4, 'MarkerSize', 12);
plot(thresh_values, MSE_or,   '-d', 'LineWidth', 4, 'MarkerSize', 12);
plot(thresh_values, MSE_os,   '-^', 'LineWidth', 4, 'MarkerSize', 12);

xlabel('$\kappa$', 'Interpreter', 'latex', ...
       'FontSize', 30, 'FontWeight','bold');
ylabel('MSE', 'FontSize', 30, 'FontWeight','bold');

set(gca,'FontSize',28,'FontWeight','bold');

legend('Graph-QMF', 'GraphBior', 'Critically-sampled','Oversampled', ...
       'FontSize', 24, 'Location','northeast');

grid on;

ax = gca;
ax.FontSize = 20;
ax.LineWidth = 1.5;