function [Ln_bpt, c] = QMFfb_filters(theta, bptG, filterlen, arange)

N = size(bptG,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Normalized Laplacian Matrices for Each Bpt graph
% Compute Normalized Laplacian Matrices for Each Bpt graph
disp('Computing normalized Laplacian matrices for each subgraph...');
Ln_bpt = zeros(N,N,theta);
for i = 1:theta
    d1 = sum(bptG(:,:,i),2);
    d1(d1 == 0) = 1; % for isolated nodes
    d1_inv = d1.^(-0.5);
    D1_inv = diag(d1_inv);
    An = D1_inv*bptG(:,:,i)*D1_inv;
    An = 0.5*(An + An');
    Ln_bpt(:,:,i) = eye(N) - An;
end

% design a low-pass kernel
% S = sprintf('%s%d%s ' , ' Computing a ',filterlen,' ^th order approximation of Meyer kernel');
% disp(S)
%g = @(x)(ideal_kernel(x)); % Ideal Wavelet Kernel
g = @(x)(meyer_kernel(x));  % Meyer Wavelet Kernel
c=sgwt_cheby_coeff(g,filterlen,filterlen+1,arange);