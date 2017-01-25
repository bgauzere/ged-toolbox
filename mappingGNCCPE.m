function [mapping mapping_time zeta] = mappingGNCCPE(G1,G2, costs, params)
    n = size(G1,1);
    m = size(G2,1);
    if (params.method == 1) % random walks sspr
        CM = RandomWalksCostMatrixLSAPE(G1,G2,costs,params.k);
        [Minit,~,~] = hungarianLSAPE(CM);
        Minit = double(Minit);
    % elseif (params.method == 2) %% Bunke
    %     C = NodeCostMatrix(G1, G2, costs);
    %     %Mapping init
    %     [phi_Minit,u,v] = hungarianLSAP(C);
    %     nplusm = length(phi_Minit);
    %     Minit=zeros(nplusm,nplusm);
    %     Minit(sub2ind([nplusm,nplusm],int32([1:nplusm])',phi_Minit+1)) = 1;%Computation of assignement matrix from phi
    % elseif (params.method == 3) %random
    %     Minit=eye(n,m);
    %     idx=randperm(n+m);
    %     Minit=Minit(idx,:);
    % elseif(params.method == 6) %just node matrix
    %     [mapping] = mappingLSAP(G1,G2,costs,params.k,6);
    %     sizeAssign = size(mapping,1);
    %     Minit = zeros(n+m,n+m);
    %     for c = 1:sizeAssign
    %         Minit(c, mapping(c)) = 1;
    %     end
    % elseif(params.method == 7) %zeros matrix
    %     Minit = eye(n+m);
    % elseif (params.method == 8) %close to minima matrix
    %     Minit = ones(n+m)./(n+m);
    end
    [map_matrix mapping_time] = gnccpe(G1, G2,costs,params.maxIter,Minit,params.d,params.debug);
    [mapping, phi_i]= find(int32(map_matrix)');
end

