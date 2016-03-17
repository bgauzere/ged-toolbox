function CM = BunkeCostMatrix(G1,G2,costs)
n = size(G1,1);
m = size(G2,1);
CM=inf(n+m,n+m);
CM(n+1:end,m+1:end) = 0;
for j=1:m
    %Neighbours of j in g2
    a2{j} = find(G2(j,:));
    a2{j}((a2{j}==j)) = [];
    edges2{j}=G2(a2{j},j);
    nodes2{j}=diag(G2(a2{j},a2{j}));
    sub_m{j}=length(a2{j});
    
    e2_mat{j} = edges2{j}(:,ones(1,n))';
    n2_mat{j} = nodes2{j}(:,ones(1,n))';
end


for i=1:n
    %Neighbours of i in g1
    a1=find(G1(i,:));
    a1((a1==i)) = [];
    edges1=G1(a1,i);
    nodes1 = diag(G1(a1,a1));
    sub_n=length(a1);
    
    e1_mat = edges1(:,ones(1, m));
    n1_mat = nodes1(:,ones(1, m));
    
    sub_costR=diag(ones(sub_n,1)*(costs.cnd + costs.ced));
    sub_costR(sub_costR == 0) = inf;

    for j=1:m
        sub_costA= ((e2_mat{j}(1:sub_n,:) ~= e1_mat(:,1:sub_m{j})).*costs.ces)  + ((n2_mat{j}(1:sub_n,:) ~= n1_mat(:,1:sub_m{j})).*costs.cns);
        sub_costI=diag(ones(sub_m{j},1)*(costs.cni + costs.cei));
        sub_costI(sub_costI == 0) = inf;
        sub_cost = [sub_costA,sub_costR; sub_costI,zeros(sub_m{j},sub_n)];
        [mapping,~,~] = hungarianLSAP(sub_cost);
        CM(i,j) = sum(sum(sub_cost((mapping'*(sub_n+sub_m{j}))+int32([1:length(mapping)])))) + ...
                  costs.cns*(G1(i,i) ~= G2(j,j));
        CM(n+j,j) = length(a2{j})*(costs.cni + costs.cei) + costs.cni;
        CM(i,m+i) = length(a1)*(costs.cnd + costs.ced) + costs.cnd;
    end
end
