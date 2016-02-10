function [ H,L ] = histoLab(nbLab, LN, W)
% [H,L] = histoLab(nbLab, LN, W )
%   compute, for each node of the graph and for each label, the number of
%   walks ending at a node with the label
%
%   nbLab : number of labels (1st label must be 1)
%   LN : label-node matrix such that LN(1,i) is the label of node i
%   W : node-node matrix representing walks

    L = zeros(nbLab,size(W,1));
    for i=1:nbLab
        L(i,:) = (LN == i);
    end
    H = L * W;
end
