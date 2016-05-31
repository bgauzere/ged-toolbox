function [mapping mapping_time zeta] = mappingGNCCP(G1,G2, costs, params)
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
    elseif(params.method == 6) %just node matrix
        [mapping] = mappingLSAP(G1,G2,costs,params.k,6);
        sizeAssign = size(mapping,1);
        Minit = zeros(n+m,n+m);
        for c = 1:sizeAssign
            Minit(c, mapping(c)) = 1;
        end
    end
    [map_matrix mapping_time zeta] = gnccp(G1, G2,costs,params.maxIter,Minit,params.d,0);
    

    %% Code pour supprimer les incohérences sur la zone epsilon
    %mapbis=zeros(size(map_matrix));
    % mapbis(find(map_matrix == repmat(max(map_matrix), ...
    %                                  size(map_matrix,1),1))) = 1;
    % map_matrix = mapbis;
    % map_matrix(n+1:end,m+1:end) = 0;
    % for (i=1:n)
    %     if (sum(mapbis(n+i,m+1:end)) > 0.1)
    %         map_matrix(n+i,m+i) = 1;
    %     end
    % end
     % map_matrix
    %int32(map_matrix)
    [mapping, phi_i]= find(int32(map_matrix)');
end

