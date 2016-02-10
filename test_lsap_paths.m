path_dir = './costpaths';
walk_dir = './costwalks';
utils_dir = './utils';
addpath(path_dir,walk_dir,utils_dir);

dataset='../data/dataset_bps.ds';

X=DatasetToAdjacency(dataset);
N = numel(X)

costs.cei = 3;
costs.ces = 1;
costs.ced = 3;
costs.cni = 3;
costs.cns = 1;
costs.cnd = 3;

ed_lsap_paths = zeros(183,183);

for i = 1:183
    for j = 1:183
        G1 = X(i).am;
        G2 = X(j).am;
        n = size(G1,1);
        idx=randperm(n);
        G1=G1(idx,idx);
        m = size(G2,1);
        idx=randperm(m);
        G2=G2(idx,idx);
        
        % First matrix is not bigger than the second one
        if(size(G2)  < size(G1))
            tmp = G1;
            G1 = G2;
            G2 = tmp;
        end
        fprintf('calcul %d,%d \n', i,j);
        params.framework =  2;
        params.k = 3;
        params.method = 2;
        [val,map] = editDistance(X(i).am,X(j).am, costs, params);
        ed_lsap_paths(i,j) = val;
    end
end
