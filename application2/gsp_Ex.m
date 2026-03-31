function [G] = gsp_Ex(N)


if nargin < 1
   N = 989; 
end

G.N = N;

load Ex_graph.mat


% Create weighted adjacency matrix
i_inds = 1:N;                  % indices of each node
j_inds = [N, 1:(N-1)];         % indices of the next node (cyclic)

% Only one direction is added, from node i to node i+1 (circular)、
G.W = A;

% Create coordinates
G.coords = xy;
G.plotting.limits = [-1, 1, -1, 1];

G.type = 'Ex';

G = gsp_graph_default_parameters(G);

end