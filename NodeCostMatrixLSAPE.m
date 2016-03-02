%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  NodeCostMatrix: compute the cost matrix between nodes
%  based on parameters:
%           A1       first  adjacency matrix 
%           A2       second  adjacency matrix
%           costs
%
% NB: The matrix is squared, having a dimnesion equals to size(A1)+size(A2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = NodeCostMatrixLSAPE(A1, A2, costs)

% Substitution parts
nodes_A1 = diag(A1);
n=length(diag(A1));

nodes_A2 = diag(A2);
m=length(diag(A2));

A_part1 = (repmat(nodes_A1,1,m) ~= repmat(nodes_A2',n,1))*costs.cns;

%Matrix association
A = [A_part1, ones(n,1).*costs.cnd ; ones(1,m).*costs.cni, 0];
return
