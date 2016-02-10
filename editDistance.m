function [edit_distance, mapping] = editDistance(G1,G2, costs, ...
                                                 params)
% Mettre le controle de la présence des params
if (params.framework == 1)
    mapping = mappingQAP(G1,G2,costs,params);
elseif (params.framework == 2)
    mapping = mappingLSAP(G1,G2,costs, params.k , params.method);
end
edit_distance = computeEditDistance(G1,G2,int32(mapping),costs.cns, ...
                                    costs.cnd,costs.ces,costs.ced);
