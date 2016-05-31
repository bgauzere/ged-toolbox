function [proj_map,  running_time,zeta] = gnccp(G1, G2, costs,maxIter, ...
                                           Mapping_init, d, debug)
zeta = 1;
if nargin < 7
    debug = 0;
    if nargin < 6
        d = 0.01;
    end
end

mapping = Mapping_init;
flag = true;
n = size(G1,1);
m = size(G2,1);
while ((zeta > -1) && flag)
    % if (zeta < 0.030001)
    %     debug =  1;
    % end
    if debug 
        fprintf('Iteration ipfp avec zeta = %f \n', zeta); 
    end

    [mapping, running_time,conv,S_zeta] = ipfp_zeta(G1, G2,costs,maxIter, ...
                                                    mapping,zeta, ...
                                                    debug);    
    
    flag_n = sum(abs(max(mapping(1:n,:)) - sum(mapping(1:n,:))))/n > 0.001 ;
    flag_m = sum(abs(max(mapping(:,1:m),[],2) - sum(mapping(:,1:m),2)))/m > 0.001 ;
    flag = flag_n * flag_m;
    zeta = zeta - d;
    if debug   
        plot(S_zeta)
        % sum(abs(max(mapping) - sum(mapping)))
        mapping
        proj_map =zeros(size(mapping));
        [~, I] = max(mapping);
        proj_map(sub2ind([n+m,n+m],I,1:n+m)) = 1;
        proj_map
        int32(mapping)
        pause;
    end
    proj_map =zeros(size(mapping));
    [~, I] = max(mapping);
    proj_map(sub2ind([n+m,n+m],I,1:n+m)) = 1;
end
