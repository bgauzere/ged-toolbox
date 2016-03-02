function [mapping, mapping_time] = mappingLSAPE(G1,G2, costs, kw,method)
    n = size(G1,1);
    m = size(G2,1);
    
    G1=double(G1);
    G2=double(G2);
    
    if (method == 1) %Original (SSPR 2014)
        % CM = RandomWalksCostMatrixLSAPE(G1,G2,costs,kw);
        CM = RandomWalksCostMatrixLSAPE(G1,G2,costs,kw);
    elseif (method  == 2) % Paths
        BoB1 = bagOfBagsOfSimpleLabeledPaths(G1,kw);
        BoB2 = bagOfBagsOfSimpleLabeledPaths(G2,kw);
        [CM] = costBagsOfBagsOfSimpleLabeledPaths(BoB1,BoB2, ...
                                                  costs.cns,costs.cnd,costs.ces,costs.ced,2,@ ...
                                                  costLabeledPaths); ...
        % Luc
        CM=LSAPtoLSAPECostMatrix(CM,n,m);
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
        0;
    elseif (method == 5) % random
        S = randi(20,n,m);
        D = randi(20,n,1);
        A  = randi(20,1,m);
        CM = [ S D ; A 0 ];
    elseif(method == 6) %just node matrix
        CM = NodeCostMatrixLSAPE(G1, G2, costs);
    end
    chrono_mapping=tic;
    [sr1,sc1,u,v] = hungarianLSAPE(CM);
    mapping_time = toc(chrono_mapping);
    %convert sr1 and sc1 to original mapping format
    mapping=LSAPEtoLSAPMapping(sr1,sc1);
end
