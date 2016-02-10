function [CM] = RandomWalksCostMatrix(G1,G2,costs,kw)

n = size(G1,1);
m = size(G2,1);
nbLab = max(max(diag(G1)),max(diag(G2)));

[C,Cie,Cej] = costLabeledWalksApprox(nbLab,kw,G1,G2,costs.cns,costs.ces,costs.cnd,costs.ced);

Ti = inf*ones(n);
Ti(logical(eye(n))) = Cie';
Tj = inf*ones(m);
Tj(logical(eye(m))) = Cej';
CM = [[reshape(C,m,n)',Ti];[Tj,zeros(m,n)]];

return
