function [F,V,L] = gft(W)

[V,L] = eig(full(W)); % Decomposition of the adjacency matrix
[~,Index] = sort(diag(L),'descend'); % Sort
V = V(:,Index);
F = inv(V);
end