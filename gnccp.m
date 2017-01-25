function [proj_map,  running_time] = gnccp(G1, G2, costs,maxIter, ...
                                           Mapping_init, d, debug)
zeta = 1;
if nargin < 7
    debug = 0;
    if nargin < 6
        d = 0.01;
    end
end
precision = 1e-8;
mapping = Mapping_init;
flag = true;
n = size(G1,1);
m = size(G2,1);

while ((zeta > -1) && flag)
    if debug 
        fprintf('Iteration ipfp avec zeta = %f \n', zeta); 
    end

    [mapping, running_time] = ipfp_zeta(G1, G2,costs,maxIter, ...
                                                    mapping,zeta, ...
                                                    0);    
    flag = (sum(sum(mapping - double(int32(mapping)))) > precision);
    zeta = zeta - d;
    if debug   
        disp(sum(sum(mapping - double(int32(mapping)))));
        proj_map =zeros(size(mapping));
        [~, I] = max(mapping');
        proj_map(sub2ind([n+m,n+m],I,1:n+m)) = 1;
        [formatted_mapping, phi_i]= find(int32(proj_map)');
        edit_distance = computeEditDistance(G1,G2,int32(formatted_mapping),costs.cns, ...
                                            costs.cnd,costs.ces, ...
                                            costs.ced);
    end
end
[phi_sub_problem,~,~] = hungarianLSAP(1-mapping); %b= argmax x^tb
phi_sub_problem+1;
proj_map=zeros(n+m,n+m);
proj_map(sub2ind([n+m,n+m],int32((1:n+m)'),phi_sub_problem+1)) = 1;%Computation of assignement matrix from phi
