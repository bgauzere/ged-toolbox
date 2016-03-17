function [mapping, mapping_time] = mappingLSAP(G1,G2, costs, kw,method)
    nbLab = max(max(diag(G1)),max(diag(G2)));
    n = size(G1,1);
    m = size(G2,1);
    
    G1=double(G1);
    G2=double(G2);
    
    if (method == 1) %Original (SSPR 2014)
        CM = RandomWalksCostMatrix(G1,G2,costs,kw);
    elseif (method  == 2) % Paths
        BoB1 = bagOfBagsOfSimpleLabeledPaths(G1,kw);
        BoB2 = bagOfBagsOfSimpleLabeledPaths(G2,kw);
        [CM] = costBagsOfBagsOfSimpleLabeledPaths(BoB1,BoB2, ...
                                                  costs.cns,costs.cnd,costs.ces,costs.ced,2,@ ...
                                                  costLabeledPaths); % Luc
    elseif (method == 3) % ipfp
        0;
        % elseif (method  == 5) Gbr 2015
        %   [CM] = computeExactGED(G1,G2, costs.cns, costs.ces, costs.cnd, costs.ced, kw);
        % elseif (method  == 6)
        %   [CM] = computeApproxKGraphsGed(G1,G2, costs.cns, costs.ces, costs.cnd, costs.ced, kw,limit);
        % elseif (method  == 7)
        %   [CM] = computeApproxKGraphsHGGed(G1,G2, costs.cns,
        %   costs.ces, costs.cnd, costs.ced, kw);
    elseif (method == 4) %Bunke/Riesen
        CM=BunkeCostMatrix(G1,G2,costs);
    elseif (method == 5) % random
        S = randi(20,n,m);
        D = diag(randi(20,n,1));
        D(D==0) = inf;
        A  = diag(randi(20,m,1));
        A(A==0) = inf;
        CM = [ S D ; A zeros(m,n) ]; 
    elseif(method == 6) %just node matrix
        CM = NodeCostMatrix(G1, G2, costs);
    end
    chrono_mapping=tic;
    [mapping,u,v] = hungarianLSAP(CM);
    mapping_time = toc(chrono_mapping);
    mapping = mapping +1;
end
