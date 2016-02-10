function [ c ] = costLabeledPaths( s1, s2, cns, cnr, ces, cer )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

    m = (numel(s1)-1)/2;
    n = (numel(s2)-1)/2;
    
    D = zeros(m+1,n+1);
    cr = cnr + cer;

    for i=2:m+1
        D(i,1) = D(i-1,1) + cr;
    end
    
    for i=2:n+1
      D(1,i) = D(1,i-1) + cr;
    end
    
    for j=2:n+1
        for i=2:m+1
            if (s1(2*(i-1)+1,1) == s2(2*(j-1)+1,1)) % same node label
                if (s1(2*(i-1),1) == s2(2*(j-1),1)) % same edge label
                    D(i,j) = D(i-1,j-1);
                else
                    D(i,j) = ces;
                end
            else % different node label
                if (s1(2*(i-1),1) == s2(2*(j-1),1)) % same edge label
                    D(i,j) = min([D(i-1,j)+cr,D(i,j-1)+cr,D(i-1,j-1)+cns]);
                else
                    D(i,j) = min([D(i-1,j)+cr,D(i,j-1)+cr,D(i-1,j-1)+cns+ces]);
                end
            end
        end
    end
    
    c = D(m+1,n+1);

end
