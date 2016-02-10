function [ c] = costBagsOfSimpleLabeledPaths(Bi, k, Bj, l, cns, cnr, ces, cer)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    
    % -------------------------------------------------------------
    % different levels
    if (k ~= l)
       
        kl = min(k,l);
        % cost of assigning non-empty bags having same maximal path length
        c = costBagsOfSimpleLabeledPathsEqual(Bi{kl}, Bj{kl}, cns, cnr, ces, cer);
        
        if kl == k
            for i=k:l-1
                c = c + (cnr+cer)*(size(Bj{i+1},2) + (size(Bj{i+1},2)-size(Bj{i},2))*i) + cnr;
            end
        else
            for i=l:k-1
                c = c + (cnr+cer)*(size(Bi{i+1},2) + (size(Bi{i+1},2)-size(Bi{i},2))*i) + cnr;
            end
        end
        
        return;
    end
    
    % -------------------------------------------------------------
    % same level (same number of nodes), that is k == l
    c = costBagsOfSimpleLabeledPathsEqual(Bi{k}, Bj{l}, cns, cnr, ces, cer);
    
end
