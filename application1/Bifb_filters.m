function [Ln_bpt, c_d, c_r] = Bifb_filters(theta, bptG, filterlen, arange)

N = size(bptG,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% design biorthogonal kernels
N_lo = floor(filterlen/2) + mod(filterlen,2);
N_hi = N_lo;
% S = sprintf('%s%d%s ' , ' Computing a ',filterlen,' ^th order approximation of Meyer kernel');
% disp(S)
[hi_d,lo_d] = biorth_kernel(N_lo,N_hi); % Nd zeros for lowpass Nr zeros for highpass
filterlen_hi = length(roots(hi_d));
filterlen_lo = length(roots(lo_d));
g0 = @(x)(polyval(lo_d,x));
g1 = @(x)(polyval(hi_d,x));
h0 = @(x)(polyval(hi_d,2 - x));
h1 = @(x)(polyval(lo_d,2 - x));

c_d{1}=sgwt_cheby_coeff(h0,filterlen_lo,filterlen_lo+1,arange);
c_d{2}=sgwt_cheby_coeff(h1,filterlen_hi,filterlen_hi+1,arange);
c_r{1}=sgwt_cheby_coeff(g0,filterlen_hi,filterlen_hi+1,arange);
c_r{2}=sgwt_cheby_coeff(g1,filterlen_lo,filterlen_lo+1,arange);