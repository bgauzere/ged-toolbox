%Your path to ged toolbox
ged_toolbox_path='/home/bgauzere/dev/ged-framework/ged-toolbox/';
path_dir = [ged_toolbox_path 'costpaths'];
walk_dir = [ged_toolbox_path 'costwalks'];
utils_dir = [ged_toolbox_path 'utils'];
addpath(path_dir,walk_dir,utils_dir,ged_toolbox_path);

%Your path to Synthetic dataset location
dataset_path='/home/bgauzere/dev/ged-framework/data/';


costs.cei = 3;
costs.ces = 1;
costs.ced = 3;
costs.cni = 3;
costs.cns = 1;
costs.cnd = 3;

sizes = [10,20,25,30,40,50,60,70,80,90,100,250,500];
N=100;
distance_synthetic_lsape = zeros(N,length(sizes));
mapping_times_lsape = zeros(length(sizes),1);

for cursize = 1:length(sizes)
    fprintf('Size : %d \n', sizes(cursize));
    for i = 1:N
        Source=DatasetToAdjacency([dataset_path, 'data-random-',num2str(sizes(cursize)), '-source.ds']);
        Target= DatasetToAdjacency([dataset_path,'data-random-',num2str(sizes(cursize)), '-target.ds']);
        
        G1 = Source(i).am;
        G2 = Target(i).am;
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
        fprintf('Computation of pair %i for size %d \n',i,sizes(cursize));
        params.framework =  3;
        params.k = 3;
        params.method = 6;
        [val,map,mapping_time] = editDistance(G1,G2, costs, params);
        mapping_times_lsape(cursize) =  mapping_times_lsape(cursize) + mapping_time;
        distance_synthetic_lsape(i,cursize) = val;
    end
end
