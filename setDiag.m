function [ M ] = setDiag( M, D )
%setDiag(M,D) set D to the diagonal of matrix M
%   
    M(logical(eye(size(M)))) = D;

end

