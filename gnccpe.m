function [proj_map,  running_time] = gnccpe(G1, G2, costs,maxIter, ...
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

    [mapping, running_time] = ipfpe_zeta(G1, G2,costs,maxIter, ...
                                                    mapping,zeta, ...
                                         0);    
    % flag_n = sum(abs(max(mapping(1:n,:)) - sum(mapping(1:n,:))))/n > 0.001 ;
    % flag_m = sum(abs(max(mapping(:,1:m),[],2) - sum(mapping(:,1:m),2)))/m > 0.001 ;
    % flag = flag_n * flag_m;
    % if(zeta ==1)
    %     disp('mapping pour zeta = 1');
    %     disp(mapping);
    %     pause;
    % end
    flag = (sum(sum(mapping - double(int32(mapping)))) > precision);
    zeta = zeta - d;
    if debug   
        % [ S_zeta S]
        disp(sum(sum(mapping - double(int32(mapping)))));
        mapping
        % figure();
        % image(mapping*100);
        % colormap('gray');
        % pause;
        % save('mapping_sortie_gnccp.mat','mapping')
        % proj_map =zeros(size(mapping));
        % [~, I] = max(mapping');
        % proj_map(sub2ind([n+m,n+m],I,1:n+m)) = 1;
        % save('proj_map.mat','proj_map')
        % [formatted_mapping, phi_i]= find(int32(proj_map)');
        % save('formatted_mapping.mat','formatted_mapping')
        % edit_distance = computeEditDistance(G1,G2,int32(formatted_mapping),costs.cns, ...
        %                                     costs.cnd,costs.ces, ...
        %                                     costs.ced);
        
        % disp(edit_distance)
        % [integer_mapping, phi_i]= find(bkp1');
        % edit_distance = computeEditDistance(G1,G2,int32(integer_mapping),costs.cns, ...
        %                                     costs.cnd,costs.ces, ...
        %                                     costs.ced);
        % disp(edit_distance)
        
        % [integer_mapping, phi_i]= find(bkp1');
    end
end
[proj_map,~,~] = hungarianLSAPE(1-mapping); %b= argmax x^tb
