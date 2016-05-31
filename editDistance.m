function [edit_distance, mapping,mapping_time zeta] = editDistance(G1,G2, costs, ...
                                                      params)
    mapping_time=0;
    %TODO : control params
    zeta = 12;
    if (params.framework == 1)
        [mapping cost mapping_time] = mappingQAP(G1,G2,costs,params);
    elseif (params.framework == 2)
        [mapping,mapping_time] = mappingLSAP(G1,G2,costs, params.k , params.method);
    elseif (params.framework == 3)
        [mapping,mapping_time] = mappingLSAPE(G1,G2,costs, params.k ...
                                              , params.method);
    elseif (params.framework == 4)
        [mapping cost mapping_time] = mappingQAPE(G1,G2,costs,params);
    elseif (params.framework == 5)
        [mapping mapping_time zeta] = mappingGNCCP(G1,G2,costs,params);

    end
       edit_distance = computeEditDistance(G1,G2,int32(mapping),costs.cns, ...
                                        costs.cnd,costs.ces,costs.ced);
    
