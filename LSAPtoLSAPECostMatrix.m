function R = LSAPtoLSAPECostMatrix(CM,n,m)
R=zeros(n+1,m+1);
R(1:n,1:m) = CM(1:n,1:m);
R(n+1,1:end-1) = diag(CM(n+1:end,1:m));
R(1:end-1,m+1) = diag(CM(1:n,m+1:end))';
end

