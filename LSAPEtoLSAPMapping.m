function mapping = LSAPEtoLSAPMapping(sr1,sc1)
sr1 = sr1+1;
sc1 = sc1 +1;
n=length(sr1);
m=length(sc1);
mapping=zeros(n+m,1);
mapping(1:n) = sr1;
unmapped = setdiff(1:m,sr1);
mapping(n+1:n+length(unmapped)) = unmapped;
mapping(n+1+length(unmapped):end) = n+1+length(unmapped):n+m;
end
