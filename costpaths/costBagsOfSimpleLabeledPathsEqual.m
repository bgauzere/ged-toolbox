function [c] = costBagsOfSimpleLabeledPathsEqual(Bik, Bjk, cns, cnr, ces, cer)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    
    nbNodesEdges = size(Bik,1);
    nbEdges = (nbNodesEdges-1)/2;
    nbNodes = nbEdges+1;
    nbi = size(Bik,2);
    nbj = size(Bjk,2);

    % same level (same number of nodes)
    
    if (nbi == 0)
        if (nbj == 0)
            c = 0;
            return;
        end
        c = (nbNodes*cnr + nbEdges*cer)*nbj;
        return;
    end

    if (nbj == 0)
        c = (nbNodes*cnr + nbEdges*cer)*nbi;
        return;
    end

    C = zeros(nbi,nbj);
    Idxi=[];
    Idxj=[];
    for l=1:nbi
        for m=1:nbj
            for n=1:nbEdges
                if (Bik(2*n,l) ~= Bjk(2*n,m))
                    C(l,m) = C(l,m) + ces;
                end
                if (Bik(2*n-1,l) ~= Bjk(2*n-1,m))
                    C(l,m) = C(l,m) + cns;
                end
            end
            if (Bik(nbNodesEdges,l) ~= Bjk(nbNodesEdges,m))
                C(l,m) = C(l,m) + cns;
            end
            if (C(l,m) == 0)
                Idxi = [Idxi,l];
                Idxj = [Idxj,m];
            end
        end
    end
    C(Idxi,:)=[];
    C(:,Idxj)=[];
    c = 0;
    if (size(C,1)==0)
        if (size(C,2) == 0) % exact matching
            return;
        end
        c = size(C,2)*(cnr*nbNodes+cer*nbEdges); % only removal/insertion
        return;
    end
    if (size(C,2)==0)
        c = size(C,1)*(cnr*nbNodes+cer*nbEdges); % only removal/insertion
        return;
    end
    % other cases, there exists at least one column and one row in C
    % apply lsap
    cri = cnr*nbNodes + cer*nbEdges;
    Cie = cri.*ones(1,size(C,1));
    Cje = cri.*ones(1,size(C,2));
    Ti = setDiag(inf*ones(size(C,1)),Cie);
    Tj = setDiag(inf*ones(size(C,2)),Cje);
    C = [[C,Ti];[Tj,zeros(size(C,2),size(C,1))]];
    [A,c] = munkres(C);

end

