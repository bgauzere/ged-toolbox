function [ BoB ] = bagOfBagsOfSimpleLabeledPaths( M, k )
    
    BoB = cell(1,size(M,1));
    for i=1:size(M,1)
        B = {};
        [B{1:k}] = bagOfSimpleLabeledPaths(M,i,k);
        BoB{i} = B;
    end
end
