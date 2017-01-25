function [mapping cost mapping_time] = mappingQAPE(G1,G2, costs, params)
    n = size(G1,1);
    m = size(G2,1);
    if (params.method == 1) % random walks sspr
        CM = RandomWalksCostMatrixLSAPE(G1,G2,costs,params.k);
        [Minit,~,~] = hungarianLSAPE(CM);
        Minit = double(Minit);
    elseif (params.method == 2) %% Bunke
        R = NodeCostMatrixLSAPE(G1, G2, costs);
        %Mapping init
        [Minit,~,~] = hungarianLSAPE(R);
    elseif (params.method == 3) %random
        Minit=eye(n+m);
        idx=randperm(n+m);
        Minit=Minit(idx,:);
    elseif(params.method == 6) %just node matrix
        CM = NodeCostMatrixLSAPE(G1, G2, costs);
        [Minit,~,~] = hungarianLSAPE(CM);

    end
    %OK here
    Minit = double(Minit);
    [map_matrix cost mapping_time] = ipfpLSAPE(G1, G2,costs,params.maxIter, ...
                                  Minit,0);
    [mapping, phi_i]= find(map_matrix');
end
