function CM = BunkeCostMatrixLSAPE(G1,G2,costs)
n = size(G1,1);
m = size(G2,1);
CM=zeros(n+1,m+1);

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
    
    sub_costR=ones(sub_n,1).*(costs.cnd + costs.ced);
    for j=1:m
        sub_costA= ((e2_mat{j}(1:sub_n,:) ~= e1_mat(:,1:sub_m{j})).*costs.ces)  + ((n2_mat{j}(1:sub_n,:) ~= n1_mat(:,1:sub_m{j})).*costs.cns);
        sub_costI=ones(1,sub_m{j}).*(costs.cni + costs.cei);
        
        sub_cost = [sub_costA,sub_costR; sub_costI,0];
        [sr1,sc1,~,~] = hungarianLSAPE(sub_cost);
        E1 = sum(sub_cost((double(sr1)'*size(sub_cost,1) + [1:length(sr1)])));
        solLstC = find(sc1+1 == size(sub_cost,1));
        if size(solLstC,2) > 0
             E1 = E1 + sum(sub_cost(end,solLstC));
        end
        
        CM(i,j) = E1 + costs.cns*(G1(i,i) ~= G2(j,j));
        CM(n+1,j) = length(a2{j})*(costs.cni + costs.cei) + costs.cni;
        CM(i,m+1) = length(a1)*(costs.cnd + costs.ced) + costs.cnd;
    end
end
