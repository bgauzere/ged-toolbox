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

function A = NodeCostMatrix(A1, A2, costs)

%max value (compatible with C++ and Matlab)
%scaled-10000 to be compatible with Jonker approximation
maxValue = 214748;

% Substitution parts
A_part1 = ones(size(A1,1), size(A2,1))*costs.cns;
%No sub costs if nodes are equals
for i = 1:size(A1,1)
    for j = 1:size(A2,1)
        if A1(i,i) == A2(j,j)
            A_part1(i,j) = 0; 
        end
    end
end

% Deletion costs
A_part2 = ones(size(A1))*maxValue;
for i = 1:size(A1,1)
    A_part2(i,i) = costs.cnd;
end

%Insertion costs
A_part3 = ones(size(A2))*maxValue;
for i = 1:size(A2,1)
    A_part3(i,i) = costs.cni;
end

% Epsilon on epsilon mapping = zeros costs
A_part4 = zeros(size(A_part1,2), size(A_part1,1));

%Matrix association
A = [A_part1, A_part2; A_part3, A_part4];
return
