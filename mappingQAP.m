function [mapping cost] = mappingQAP(G1,G2, costs, params)
    n = size(G1,1);
    m = size(G2,1);
    if (params.method == 1) % random walks sspr
        [mapping] = mappingLSAP(G1,G2,costs,params.k,1);
        sizeAssign = size(mapping,1);
        Minit = zeros(n+m,n+m);
        for c = 1:sizeAssign
            Minit(c, mapping(c)) = 1;
        end
    elseif (params.method == 2) %% Bunke
        C = NodeCostMatrix(G1, G2, costs);
        %Mapping init
        [phi_Minit,u,v] = hungarianLSAP(C);
        nplusm = length(phi_Minit);
        Minit=zeros(nplusm,nplusm);
        Minit(sub2ind([nplusm,nplusm],int32([1:nplusm])',phi_Minit+1)) = 1;%Computation of assignement matrix from phi
    elseif (params.method == 3) %random
        Minit=eye(n+m);
        idx=randperm(n+m);
        Minit=Minit(idx,:);
    end
    [map_matrix cost] = ipfp(G1, G2,costs,params.maxIter,Minit,0);
    [mapping, phi_i]= find(map_matrix');
end

