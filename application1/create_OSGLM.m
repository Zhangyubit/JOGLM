%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Copyright (c) 2016 Akie Sakiyama
%
%  Description:
%  This function constructs an oversampled bipartite graph. 
%  
%  Inputs: 
%  A - original graph adjacency matrix
%  Colorednodes - set of colored nodes
%  Fnum - number of coloring
%
%  Outputs: 
%  bptGs - oversampled graph adjacency matrix 
%  Colorednodes - updated set of nodes (in the oversampled graph) having 
%  the same colors
%  connect_bpt - index of nodes which is duplicated in the oversampled 
%  graph
%
%  Example is given in the paper:
%% A. Sakiyama and Y. Tanaka, "Oversampled Graph Laplacian Matrix for Graph
%% Filter Banks", IEEE Transactions on Signal Processing,
%  url: http://ieeexplore.ieee.org/xpl/login.jsp?tp=&arnumber=6937182, Dec,
%  2014.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bptGs, Colorednodes, connect_bpt] = create_OSGLM(A,Colorednodes,Fnum)

N = length(A);

Colorednodes{5}=[];

for i= 1: Fnum-1;
     Colorednodes{5} = union(Colorednodes{5},Colorednodes{i});  
end
Colorednodes{6} = Colorednodes{Fnum};

A2=zeros(N);
A2(Colorednodes{5},Colorednodes{6})=A(Colorednodes{5},Colorednodes{6});
A2=(A2+A2.');
A3=A-A2;
k=sum(A3,2);
connect_bpt=find(k);
bptG1=A2;
bptG2=A3;
bptGs = cell(1);

A5=bptG2+speye(N); % add vertical edges
A5=A5(connect_bpt,:);
A6=fliplr(blkdiag(fliplr(A5.'),fliplr(A5)));
A6(1:N,1:N)=bptG1;
bptGs{1}=A6;

