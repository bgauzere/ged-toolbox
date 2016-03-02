%Your path to ged toolbox
ged_toolbox_path='/home/bgauzere/dev/ged-framework/ged-toolbox/';
path_dir = [ged_toolbox_path 'costpaths'];
walk_dir = [ged_toolbox_path 'costwalks'];
utils_dir = [ged_toolbox_path 'utils'];
addpath(path_dir,walk_dir,utils_dir,ged_toolbox_path);

%Your path to Synthetic dataset location
dataset_path='/home/bgauzere/dev/ged-framework/data/';
dataset=[dataset_path 'dataset_bps.ds'];

X=DatasetToAdjacency(dataset);
N = numel(X)

costs.cei = 3;
costs.ces = 1;
costs.ced = 3;
costs.cni = 3;
costs.cns = 1;
costs.cnd = 3;

ed_qap_random = zeros(N,N);

for i = 1:N
    for j = 1:N
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
        fprintf('Computation of  ged(%d,%d) \n', i,j);
       
        params.framework = 1;
        params.maxIter = 10;
        params.k = 5;
        params.method = 3;

        [val,map] = editDistance(X(i).am,X(j).am, costs, params);
        ed_qap_random(i,j) = val;
        
    end
end
