function [ C ] = costBagsOfBagsOfSimpleLabeledPaths(Bi, Bj , cns, cnr, ces, cer, t, costPathsFunc)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    nbNi = size(Bi,2);
    nbNj = size(Bj,2);
    C = zeros(nbNi,nbNj);
    nbK = size(Bi{1},2);
    
    % ------------------------------------------------------------------
    if t==0 % single level case
        
        for i=1:nbNi
            nbKi = nbK;
            while (size(Bi{i}{nbKi},2) == 0)
                nbKi = nbKi-1;
            end
            for j=1:nbNj
                nbKj = nbK;
                while (size(Bj{j}{nbKj},2) == 0)
                    nbKj = nbKj-1;
                end
                C(i,j) = costBagsOfSimpleLabeledPaths(Bi{i}, nbKi, Bj{j}, nbKj, cns, cnr, ces, cer);
            end
        end
        
        Cie = zeros(1,size(C,1));
        Cje = zeros(1,size(C,2));
        for i=1:nbNi
            nbK = size(Bi{i},2);
            for l=1:nbK % for each pair of bags of (k-1)-paths
                nbE = (size(Bi{i}{l},1)-1)/2;
                Cie(1,i) = Cie(1,i) + (cnr * (nbE+1) + cer * nbE) * size(Bi{i}{l},2);
            end
        end
        for j=1:nbNj
            nbK = size(Bj{j},2);
            for l=1:nbK % for each pair of bags of (k-1)-paths
                nbE = (size(Bj{j}{l},1)-1)/2;
                Cje(1,j) = Cje(1,j) + (cnr * (nbE+1) + cer * nbE) * size(Bj{j}{l},2);
            end
        end
        Ti = setDiag(inf*ones(size(C,1)),Cie);
        Tj = setDiag(inf*ones(size(C,2)),Cje);
        C = [[C,Ti];[Tj,zeros(size(C,2),size(C,1))]];
        return;
    end
    
    % ------------------------------------------------------------------
    if t==1 % multilevels case
    
        for i=1:nbNi
            nbK = size(Bi{i},2);
            for j=1:nbNj
                Cij = zeros(1,nbK);
                for l=1:nbK % for each pair of bags of (k-1)-paths
                    Cij(1,l) = costBagsOfSimpleLabeledPaths(Bi{i}, l, Bj{j}, l, cns, cnr, ces, cer);
                end
                C(i,j) = sum(Cij); % TODO better cost
            end
        end

        Cie = zeros(1,size(C,1));
        Cje = zeros(1,size(C,2));
        for i=1:nbNi
            nbK = size(Bi{i},2);
            for l=1:nbK % for each pair of bags of (k-1)-paths
                nbE = (size(Bi{i}{l},1)-1)/2;
                Cie(1,i) = Cie(1,i) + (cnr * (nbE+1) + cer * nbE) * size(Bi{i}{l},2);
            end
        end
        for j=1:nbNj
            nbK = size(Bj{j},2);
            for l=1:nbK % for each pair of bags of (k-1)-paths
                nbE = (size(Bj{j}{l},1)-1)/2;
                Cje(1,j) = Cje(1,j) + (cnr * (nbE+1) + cer * nbE) * size(Bj{j}{l},2);
            end
        end
        Ti = setDiag(inf*ones(size(C,1)),Cie);
        Tj = setDiag(inf*ones(size(C,2)),Cje);
        C = [[C,Ti];[Tj,zeros(size(C,2),size(C,1))]];
        return;
    end
    
    % ------------------------------------------------------------------
    if t==2  % bags of all k-paths with k varying (Luc)
       
        Cie = inf*ones(size(C,1));
        Cje = inf*ones(size(C,2));
        for i=1:nbNi
            for j=1:nbNj
                [C(i,j),Cie(i,i),Cje(j,j)] = costBagsOfSimplePaths(Bi{i},Bj{j},costPathsFunc,cns,cnr,ces,cer);
            end
        end
        C = [[C,Cie];[Cje,zeros(size(C,2),size(C,1))]];
        return;
    end
    
    % ------------------------------------------------------------------
    if t==3  % using non-recovering (Benoit)
       
        
        
    end
    
end
