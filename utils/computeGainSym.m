function [gain, reg, edit_distance_sym] =  computeGainSym(edit_distance, k_max)
    edit_distance_sym = zeros (size(edit_distance));
    gain=zeros(k_max,k_max);
    for i = 1:k_max
        edit_distance_sym(:,:,i) = min(edit_distance(:,:,i), ...
                                       edit_distance(:,:,i)');
        
    for i = 1:k_max
        for j = 1:k_max
            gain(i,j) = sum(sum((edit_distance_sym(:,:,i) - edit_distance_sym(:,:,j)) > 0));
        end
    end
    
    reg=zeros(k_max,k_max);
    for i = 1:k_max
        for j = 1:k_max
            reg(i,j) = sum(sum((edit_distance_sym(:,:,i) - edit_distance_sym(:,:,j)) < 0));
        end
    end
end
