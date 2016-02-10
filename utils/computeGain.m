function [gain, reg] =  computeGain(edit_distance, k_max)
    gain=zeros(k_max,k_max);
    for i = 1:k_max
        for j = 1:k_max
            gain(i,j) = sum(sum((edit_distance(:,:,i) - edit_distance(:,:,j)) > 0));
        end
    end
    
    reg=zeros(k_max,k_max);
    for i = 1:k_max
        for j = 1:k_max
            reg(i,j) = sum(sum((edit_distance(:,:,i) - edit_distance(:,:,j)) < 0));
        end
    end
end
