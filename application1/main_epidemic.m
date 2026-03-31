%% Joint Time-Vertex Oversampled Graph Filter Bank for Graph Signal Reconstruction
close all; clear; clc;

%% Load graph data
gsp_start;

G = gsp_minnesota();

figure()
gsp_plot_graph(G);

A = G.W;
L = compute_normalized_laplacian(A, 'sym');
[GFT,~,~] = gft(L);

%% Epidemic diffusion simulation on graph
Tim = [20 20 Inf Inf]; 
cp  = [0.004 0.08 0.004 0.08];
X = cell(1, 4);

for ii = 1:4
    contagion_prob    = cp(ii);            % Infection probability per day
    infection_length  = 6;                 % Infection duration (days)
    params.maxTime    = round(16*365/12);  % Simulation horizon (days)
    params.population = 70;                % Population per node (airport)
    params.model      = 'SEIRS';           % Epidemic model ('SEIRS', 'SIR', 'SIRS', 'SI')
    params.immunity   = Tim(ii);           % Immunity period
    params.latency    = 2;                 % Latent period
    rng(42);                              % Fix random seed for reproducibility
    patient_zero      = randi(size(A,1),3,1); % Randomly select 3 initial infected nodes

    % Run epidemic model
    [S, I, R, E] = compartmental_epidemics(sparse(A), contagion_prob, infection_length, patient_zero, params);
    fprintf('Infection simulation completed!\n');

    % Normalize infection signal
    I   = I - mean(I(:));
    X{ii} = I;
end

%% Visualization: Graph signals at initial and final time instances
f_1 = X{1}(:,1);
TT = size(X{1},2);
f_T = X{1}(:,TT);

param.vertex_size = 200;

figure('Units', 'pixels', 'Position', [100, 100, 1024, 768]);
gsp_plot_signal(G,f_1,param);
cbar = colorbar;
cbar.FontSize = 28;
cbar.FontWeight = 'bold';

figure('Units', 'pixels', 'Position', [100, 100, 1024, 768]);
gsp_plot_signal(G,f_T,param);
cbar = colorbar;
cbar.FontSize = 28;
cbar.FontWeight = 'bold';

% Graph Fourier Transform (GFT) visualization
hatf_1 = GFT * f_1;
plot_graph_signal_3D(G, abs(hatf_1));

hatf_T = GFT * f_T;
plot_graph_signal_3D(G, abs(hatf_T));

%% Visualization of oversampled bipartite graphs
plot_oversampled_bipartite_graph(G);

GT = gsp_ring(size(X{1},2));
plot_oversampled_bipartite_graph(GT);

%% Reconstruction and performance evaluation

% Initialize containers
hatX_os_all = cell(1,4);    X_rec_os_all = cell(1,4);
hatX_cs_all = cell(1,4);    X_rec_cs_all = cell(1,4);
X_rec_bi_all = cell(1,4);   X_rec_QMF_all = cell(1,4);

MSE_os_all = zeros(1,4); NMSE_os_all = zeros(1,4); SNR_os_all = zeros(1,4);
MSE_cs_all = zeros(1,4); NMSE_cs_all = zeros(1,4); SNR_cs_all = zeros(1,4);
MSE_bi_all = zeros(1,4); NMSE_bi_all = zeros(1,4); SNR_bi_all = zeros(1,4);
MSE_QMF_all = zeros(1,4); NMSE_QMF_all = zeros(1,4); SNR_QMF_all = zeros(1,4);

for ii = 1:4
    % Oversampled graph filter bank reconstruction
    [hatX_os, X_rec_os, MSE_os, NMSE_os, SNR_os] = jtv_ogfb_reconstruction(G, X{ii});
    
    % Critically-sampled graph filter bank reconstruction
    [hatX_cs, X_rec_cs, MSE_cs, NMSE_cs, SNR_cs] = jtv_gfb_reconstruction(G, X{ii});
    
    % Biorthogonal graph filter bank reconstruction
    [X_rec_bi, MSE_bi, NMSE_bi, SNR_bi] = jtv_Bifb_reconstruction(G, X{ii});
    
    % QMF graph filter bank reconstruction
    [X_rec_QMF, MSE_QMF, NMSE_QMF, SNR_QMF] = jtv_QMFfb_reconstruction(G, X{ii});

    % Store results
    hatX_os_all{ii}   = hatX_os;
    X_rec_os_all{ii}  = X_rec_os;
    MSE_os_all(ii)    = MSE_os;    
    NMSE_os_all(ii)   = NMSE_os;    
    SNR_os_all(ii)    = SNR_os;

    hatX_cs_all{ii}   = hatX_cs;
    X_rec_cs_all{ii}  = X_rec_cs;
    MSE_cs_all(ii)    = MSE_cs;    
    NMSE_cs_all(ii)   = NMSE_cs;    
    SNR_cs_all(ii)    = SNR_cs;

    X_rec_bi_all{ii}  = X_rec_bi;
    MSE_bi_all(ii)    = MSE_bi;    
    NMSE_bi_all(ii)   = NMSE_bi;    
    SNR_bi_all(ii)    = SNR_bi;

    X_rec_QMF_all{ii} = X_rec_QMF;
    MSE_QMF_all(ii)   = MSE_QMF;   
    NMSE_QMF_all(ii)  = NMSE_QMF;  
    SNR_QMF_all(ii)   = SNR_QMF;
end

%% Spectral visualization of oversampled representations
for ii = 1:length(hatX_os_all)
    figure('Units', 'pixels', 'Position', [100, 100, 1024, 768]);
    
    hatX_disp = log1p(abs(hatX_os_all{ii}));   % Logarithmic enhancement
    hatX_disp = hatX_disp ./ max(hatX_disp(:)); % Normalization
    
    imagesc(hatX_disp);
    colormap jet; colorbar;
    
    xlabel('Oversampled Time Frequency Index','FontSize',30,'FontWeight','bold');
    ylabel('Oversampled Vertex Frequency Index','FontSize',30,'FontWeight','bold');
    set(gca,'FontSize',28,'FontWeight','bold');
end

%% Performance comparison (MSE)
signal_idx = 1:4;

figure; hold on;
plot(signal_idx, MSE_os_all,  'o--','LineWidth',2,'MarkerSize',7);
plot(signal_idx, MSE_cs_all,  's-.','LineWidth',2,'MarkerSize',7);
plot(signal_idx, MSE_bi_all,  'd-.','LineWidth',2,'MarkerSize',7);
plot(signal_idx, MSE_QMF_all, 'x-','LineWidth',2,'MarkerSize',7);

set(gca,'YScale','log'); 
grid on;

xlabel('Signal'); 
ylabel('MSE (log scale)');
legend('Oversampled GFB','Critically-sampled GFB','Biorthogonal GFB','QMF GFB','Location','best');

xticks(1:4);
xticklabels({'$\mathbf{X}_1$','$\mathbf{X}_2$','$\mathbf{X}_3$','$\mathbf{X}_4$'});
set(gca,'TickLabelInterpreter','latex');

%% Tabulated results
signal_idx = (1:4)';

T = table(signal_idx(:), ...
          MSE_os_all(:), ...
          MSE_cs_all(:), ...
          MSE_bi_all(:), ...
          MSE_QMF_all(:), ...
          'VariableNames', {'Signal', ...
                            'Oversampled_GFB', ...
                            'CriticallySampled_GFB', ...
                            'Biorthogonal_GFB', ...
                            'QMF_GFB'});

disp(T);

% Display table in figure
f = figure('Units','pixels','Position',[100 100 800 200]);
uitable(f,'Data',T{:,:}, ...
           'ColumnName',T.Properties.VariableNames, ...
           'RowName',[], ...
           'FontSize',12, ...
           'FontWeight','bold', ...
           'Units','normalized', ...
           'Position',[0 0 1 1]);