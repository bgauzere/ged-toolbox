function [ C, Cie, Cej ] = costLabeledWalksApprox(nbLab, kw, Wa, Wb, cns, ces, cnd, ced)

    Wa = double(Wa);
    Wb = double(Wb);

    % direct product
    [Wx] = labeledKron(Wa,Wb);
    
    % histogram of similar walks
    ILx = diag(Wx)';
    Wtx = (Wx-diag(diag(Wx)))~=0;
    [Hx,Lx] = histoLab(nbLab,ILx,Wtx^kw);

    ILa = diag(Wa)';
    Wta = (Wa-diag(diag(Wa)))~=0;
    [Ha,La] = histoLab(nbLab,ILa,Wta^kw);
    
    ILb = diag(Wb)';
    Wtb = (Wb-diag(diag(Wb)))~=0;
    [Hb,Lb] = histoLab(nbLab,ILb,Wtb^kw);

    Hx = floor(sqrt(Hx));
    Hi = zeros(size(Hx));
    k = 0;
    for i=1:size(Ha,2)
        for j=1:size(Hb,2)
            Hi(:,k+j) = Ha(:,i)-min(min(Ha(:,i),Hb(:,j)),Hx(:,k+j));
        end
        k = k + size(Hb,2);
    end

    Hj = Hx;
    k = 0;
    for i=1:size(Ha,2)
        for j=1:size(Hb,2)
            Hj(:,k+j) = Hb(:,j)-min(min(Ha(:,i),Hb(:,j)),Hx(:,k+j));
        end
        k = k + size(Hb,2);
    end

    S = sum(min(Hi,Hj));
    Ri = sum(Hi)-S;
    Rj = sum(Hj)-S;
    ILxt = (ILx>0);

    C = ((kw-ILxt)*cns+kw*ces).*S + ((1+kw-ILxt)*cns+kw*ces).*min(Ri,Rj) + ((1+kw-ILxt)*cnd+kw*ced).*abs(Ri-Rj);
    Cie = ((kw+1)*cnd+kw*ced).*sum(Ha);
    Cej = ((kw+1)*cnd+kw*ced).*sum(Hb);

end