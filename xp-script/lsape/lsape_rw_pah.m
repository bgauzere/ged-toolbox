%Your path to ged toolbox
ged_toolbox_path='/home/bgauzere/dev/ged-framework/ged-toolbox/';
path_dir = [ged_toolbox_path 'costpaths'];
walk_dir = [ged_toolbox_path 'costwalks'];
utils_dir = [ged_toolbox_path 'utils'];
addpath(path_dir,walk_dir,utils_dir,ged_toolbox_path);
%Your path to GREYC dataset location
dataset_path='/home/bgauzere/dev/ged-framework/data/';
dataset=[dataset_path 'dataset_pah.ds'];

X=DatasetToAdjacency(dataset);
N = numel(X)

costs.cei = 3;
costs.ces = 1;
costs.ced = 3;
costs.cni = 3;
costs.cns = 1;
costs.cnd = 3;

mapping_time_lsape_pah = 0;

ed_lsape_rw_pah = zeros(N,N);

for i = 1:N
    for j = 1:N
        %matrix permutations
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
        fprintf('Computation of ged(%d,%d) \n',i,j);
        params.framework =  3;
        params.k = 3;
        params.method = 1;
        [val,map,mapping_time] = ...
            editDistance(X(i).am,X(j).am, costs, params);

       mapping_time_lsape_pah =  mapping_time_lsape_pah + mapping_time;
        
       ed_lsape_rw_pah(i,j) = val;
    end
end
