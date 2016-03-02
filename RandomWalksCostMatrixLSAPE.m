function [CM] = RandomWalksCostMatrixLSAPE(G1,G2,costs,kw)

n = size(G1,1);
m = size(G2,1);
nbLab = max(max(diag(G1)),max(diag(G2)));

[C,Cie,Cej] = costLabeledWalksApprox(nbLab,kw,G1,G2,costs.cns, ...
                                     costs.ces,costs.cnd,costs.ced);

CM = [ reshape(C,m,n)'  Cie'; Cej 0];
return
