function [edit_distance, mapping,mapping_time] = editDistance(G1,G2, costs, ...
                                                      params)
    mapping_time=0;
    %TODO : control params
    if (params.framework == 1)
        mapping = mappingQAP(G1,G2,costs,params);
    elseif (params.framework == 2)
        [mapping,mapping_time] = mappingLSAP(G1,G2,costs, params.k , params.method);
    elseif (params.framework == 3)
        [mapping,mapping_time] = mappingLSAPE(G1,G2,costs, params.k , params.method);
    end
    edit_distance = computeEditDistance(G1,G2,int32(mapping),costs.cns, ...
                                        costs.cnd,costs.ces,costs.ced);
    
