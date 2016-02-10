function [ c, cri1, cri2 ] = costBagsOfSimplePaths( B1, B2, costPathsFunc, cns, cnr, ces, cer )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

    n = 0;
    for i=1:numel(B1)
        n = n + size(B1{i},2);
    end
    m = 0;
    for i=1:numel(B2)
        m = m + size(B2{i},2);
    end
    
    C = zeros(n,m);
    Cie = inf*ones(n);
    Cje = inf*ones(m);
    
    % ---------------------------------
    pi = 0;
    pj = 0;
    
    for i=1:numel(B1) % for each path in B1
        B1i = B1{i};
        nbEdges = (size(B1i,1)-1)/2;
        cri = (cnr + cer) * nbEdges;
        for k=1:size(B1i,2) % for each path in bag B1{i} compare to all paths of B2
            pi = pi + 1;
            Cie(pi,pi) = cri;
            for j=1:numel(B2) % for each path in B2
                B2j = B2{j};
                for l=1:size(B2j,2)
                    pj = pj + 1;
                    C(pi,pj) = costPathsFunc(B1i(:,k),B2j(:,l),cns,cnr,ces,cer);                    
                end
            end
            pj = 0;
        end
    end
    % ---------------------------------
    pj = 0;
    for j=1:numel(B2) % for each path in B2
        B2j = B2{j};
        nbEdges = (size(B2j,1)-1)/2;
        cri = (cnr + cer) * nbEdges;
        for l=1:size(B2j,2)
            pj = pj + 1;
            Cje(pj,pj) = cri;                    
        end
    end
    
    % ---------------------------------
    C = [[C,Cie];[Cje,zeros(size(C,2),size(C,1))]];

    [mapping,u,v]  = hungarianLSAP(C);
    mapping = mapping +1;
    sizeAssign = size(mapping,1);
    M = zeros(n+m,n+m);
    for i = 1:sizeAssign
        M(i, mapping(i)) = 1;
    end
    c = linearCost(C,M);
    
    % ---------------------------------
    if nargout == 3
       cri1 = cnr + sum(diag(Cie));
       cri2 = cnr + sum(diag(Cje));
    end
    % ---------------------------------
    b1 = B1{1};
    b2 = B2{1};
    if b1(1,1) ~= b2(1,1)
        c = c + cns;
    end
end

