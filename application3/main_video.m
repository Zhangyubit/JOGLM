%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Copyright (c) 2025 Yu Zhang
%
%  Description:
%  This script performs video denoising using a joint time-vertex
%  oversampled graph Laplacian framework.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear;
gsp_start

%% ================= Parameter Settings =================
sigma = 10;                % Noise standard deviation
filterlen = 12;            % Filter length
norm_type = 'sym';         % Type of Laplacian normalization

% (Optional) Convert video format to a MATLAB-readable format using ffmpeg
% system('/usr/local/bin/ffmpeg -i input.avi -c:v mjpeg -q:v 3 -an output.avi');

videoFile = 'person_fixed.avi';
v = VideoReader(videoFile);
disp(v);

fps = v.FrameRate;         % Frame rate (typically 30 fps)
duration = 3;              % Use the first 3 seconds
numFrames = floor(fps * duration);

N = 35;                    % Resize each frame to N × N
frames = zeros(N*N, numFrames);  % Each column stores one frame (vectorized)

%% ================= Step 1: Video Loading & Preprocessing =================
for k = 1:numFrames
    frameRGB = readFrame(v);                     % Read RGB frame
    frameGray = rgb2gray(frameRGB);              % Convert to grayscale
    frameResized = imresize(frameGray, [N N]);   % Resize frame
    frames(:,k) = double(frameResized(:));       % Vectorize and store
end

% Original signal matrix (each column is a frame)
signal_matrix_original = frames;
f_or = signal_matrix_original(:);                % Vectorized ground truth

% Add Gaussian noise
rng(42);                                        % Fix random seed
signal_matrix_noisy = signal_matrix_original + sigma * randn(size(signal_matrix_original));
f_noise = signal_matrix_noisy(:);

%% ================= Step 2: Graph Construction =================
G = gsp_2dgrid(N);    % Spatial graph: each frame treated as an N×N grid graph

%% ================= Step 3: Denoising =================
T = numFrames;        % Number of time samples (frames)
thresh = 4;           % Threshold parameter

%% ---- Oversampled Graph Filter Bank ----
tic;
X_rec_os = jtv_ogfb_denoise(G, signal_matrix_noisy, ...
    filterlen, norm_type, sigma, thresh);
time_os = toc;

%% ---- Critically Sampled Graph Filter Bank ----
tic;
X_rec_or = jtv_gfb_denoise(G, signal_matrix_noisy, ...
    filterlen, norm_type, sigma, thresh);
time_or = toc;

fprintf('Running time (Oversampled): %.4f seconds\n', time_os);
fprintf('Running time (Critically sampled): %.4f seconds\n', time_or);

% Dimension check
n = N * N;
if size(X_rec_os,1) ~= n || size(X_rec_os,2) ~= T
    error('X_rec_os has incorrect size, expected [N^2, T]');
end
if size(X_rec_or,1) ~= n || size(X_rec_or,2) ~= T
    error('X_rec_or has incorrect size, expected [N^2, T]');
end

%% ================= Step 4: Vectorization =================
f_hat_os = X_rec_os(:);
f_hat_or = X_rec_or(:);

%% ================= Step 5: Error Metrics =================
MSE_noisy = mean((f_or - f_noise).^2);
SNR_noisy = 10*log10(mean(f_or.^2)/MSE_noisy);

MSE_os = mean((f_or - f_hat_os).^2);
SNR_os = 10*log10(mean(f_or.^2)/MSE_os);

MSE_or = mean((f_or - f_hat_or).^2);
SNR_or = 10*log10(mean(f_or.^2)/MSE_or);

fprintf('sigma value = %s\n', rats(sigma));
fprintf('MSE of noisy signal = %.4f, SNR = %.2f dB\n', MSE_noisy, SNR_noisy);
fprintf('MSE of denoised signal by Oversampling = %.4f, SNR = %.2f dB\n', MSE_os, SNR_os);
fprintf('MSE of denoised signal by Critically Sampling = %.4f, SNR = %.2f dB\n', MSE_or, SNR_or);

% Empirically optimal thresholds: ~4 and ~1.8

%% ================= Step 6: Visualization (Frame Comparison) =================
frameIdx = [1, 10, 20, 30, 40, 50, 60, 70];  
numFrames = length(frameIdx);

figure('Position', [100 100 1800 900]);
t = tiledlayout(4,numFrames,'TileSpacing','compact','Padding','compact');

for i = 1:numFrames
    j = frameIdx(i);

    % Original frame
    nexttile(i);
    imshow(reshape(signal_matrix_original(:,j), [N,N]), []);
    title(sprintf('Original (Frame %d)', j), 'FontSize', 14, 'FontWeight','bold');

    % Noisy frame
    nexttile(i+numFrames);
    imshow(reshape(signal_matrix_noisy(:,j), [N,N]), []);
    title(sprintf('Noisy (Frame %d)', j), 'FontSize', 14, 'FontWeight','bold');

    % Critically sampled reconstruction
    nexttile(i+2*numFrames);
    imshow(reshape(X_rec_or(:,j), [N,N]), []);
    title(sprintf('Critically sampled (Frame %d)', j), 'FontSize', 14, 'FontWeight','bold');

    % Oversampled reconstruction
    nexttile(i+3*numFrames);
    imshow(reshape(X_rec_os(:,j), [N,N]), []);
    title(sprintf('Oversampled (Frame %d)', j), 'FontSize', 14, 'FontWeight','bold');
end

%% ================= SNR vs Threshold Analysis =================
thresh_values = 0:0.2:4;
num_thresh = length(thresh_values);

SNR_os_curve = zeros(num_thresh,1);
SNR_or_curve = zeros(num_thresh,1);

% Noisy SNR is constant
SNR_noisy_curve = SNR_noisy * ones(num_thresh,1);

for i = 1:num_thresh
    
    thresh = thresh_values(i);
    
    % Oversampled
    X_rec_os = jtv_ogfb_denoise(G, signal_matrix_noisy, ...
        filterlen, norm_type, sigma, thresh);
    f_hat_os = X_rec_os(:);
    MSE_os = mean((f_or - f_hat_os).^2);
    SNR_os_curve(i) = 10*log10(mean(f_or.^2)/MSE_os);
    
    % Critically sampled
    X_rec_or = jtv_gfb_denoise(G, signal_matrix_noisy, ...
        filterlen, norm_type, sigma, thresh);
    f_hat_or = X_rec_or(:);
    MSE_or = mean((f_or - f_hat_or).^2);
    SNR_or_curve(i) = 10*log10(mean(f_or.^2)/MSE_or);
    
end

%% ================= Plot SNR vs Threshold =================
figure('Position', [200, 200, 900, 500]); 
hold on;

% Color definitions
color_noisy = [0.2 0.2 0.2];   % Dark gray
color_or    = [0 0.2 0.6];     % Dark blue
color_os    = [0.6 0 0];       % Dark red

plot(thresh_values, SNR_noisy_curve, '--', ...
    'Color', color_noisy, 'LineWidth', 3);

plot(thresh_values, SNR_or_curve, '-o', ...
    'Color', color_or, 'LineWidth', 3, ...
    'MarkerSize', 8, 'MarkerFaceColor', color_or);

plot(thresh_values, SNR_os_curve, '-s', ...
    'Color', color_os, 'LineWidth', 3, ...
    'MarkerSize', 8, 'MarkerFaceColor', color_os);

xlabel('$\kappa$', 'Interpreter', 'latex', ...
    'FontSize', 20, 'FontWeight','bold');

ylabel('SNR (dB)', ...
    'FontSize', 20, 'FontWeight','bold');

set(gca,'FontSize',20,'FontWeight','bold');

legend('Noisy signal','Critically-sampled','Oversampled', ...
    'FontSize', 20, 'Location','northwest');

grid on;

ax = gca;
ax.LineWidth = 1.5;
box on;