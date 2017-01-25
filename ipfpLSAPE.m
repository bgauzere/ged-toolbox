function [mapping, current_ged, running_time] = ipfpLSAPE(G1, G2, costs, maxIter, Mapping_init, debug)
    chrono=tic;
    if nargin < 6
        debug = 0;
    end
    n=size(G1,1);
    m=size(G2,1);
    CostMatrix = NodeCostMatrixLSAPE(G1, G2, costs);
    Xk = Mapping_init; 
    [phi_i, i] = find(Xk');
    [XkD] = quadraticTermLSAPE(G1, G2, costs.cei, costs.ced, costs.ces, Xk, ...
                             [i, phi_i]);
    Lterm = linearCost(CostMatrix, Xk);
    S(1) = sum(sum(XkD.*Xk)) + Lterm;
    k=1;
    Xkm1 = 2.*Xk; % For the first loop
    while (  (k < maxIter)  && (norm(Xk - Xkm1,1) > 0.0001));
        [phi_i, i] = find(Xk');
        XtD = quadraticTermLSAPE(G1, G2, costs.cei, costs.ced, costs.ces, Xk, ...
                            [i, phi_i]);
        linear_sub_problem = 2*XtD + CostMatrix;
        % linear_sub_problem(linear_sub_problem > 100000) = -1; % passage
        %                                                       % des inf à -1
        % R = LSAPtoLSAPECostMatrix(linear_sub_problem,n,m);
        [bkp1,~,~] = hungarianLSAPE(linear_sub_problem);
        bkp1 = double(bkp1);
        R(k+1) = linearCost(linear_sub_problem,bkp1);
        Lterm_new  = linearCost(CostMatrix, bkp1);
        
        [phi_i, i] = find(bkp1');
        Xkp1tD = quadraticTermLSAPE(G1, G2, costs.cei, costs.ced, costs.ces, bkp1, ...
                               [i, phi_i]);
        quadratic_part = sum(sum( Xkp1tD.*bkp1));
        S(k+1) =  quadratic_part + Lterm_new;
        current_ged = quadratic_part/2 + Lterm_new;

        alpha = R(k+1) - 2*S(k) + Lterm;
        % alpha=sum(sum((2*XtD).*(bkp1-Xk))) + linearCost(CostMatrix,(bkp1- ...
        %                                                   Xk));
        
        beta_ = S(k+1) + S(k) - R(k+1) - Lterm;
        % [phi_i, i] = find((bkp1 - Xk)');
        % bkp1mXktD = quadraticTerm(G1, G2, costs.cei, costs.ced, costs.ces, (bkp1 - Xk), ...
        %                           [i, phi_i]);
        % beta_ = sum(sum( bkp1mXktD.*(bkp1 - Xk)));
        
        t0 = 0;
        if (beta_ > 0.00001)
            t0 = -alpha / (2* beta_); %% Compute t0 here to avoid division
            %% by 0 
        end 
        if(debug)
           fprintf('tracosts.ce : %f \n', trace(linear_sub_problem'*(bkp1-Xk)));
           fprintf('R(k+1) : %f \n',R(k+1));
           fprintf('Lterm : %f \n',Lterm);
           fprintf('S(k) : %f \n',S(k));
           fprintf('S(k+1) %f : \n', S(k+1));
            % ged_real = computeEditDistance(G1, G2, bkp1, costs.ces,
            % costs.ced, costs.ces, costs.ced); //Refaire le mapping
            % plutot que bkp1 pour que ça marche
            % fprintf('real ged : %f \n',ged_real); 
            sum(sum( Xkp1tD.*bkp1))
            linearCost(CostMatrix, bkp1)
           fprintf('my current ged %f : \n', sum(sum( Xkp1tD.*bkp1))/2 + linearCost(CostMatrix, bkp1));
            %fprintf('2*X^tD(bkp1 - xk) + c^t(bkp1 - xk) : %f \n',sum(sum((2*XtD).*(bkp1-Xk))) ...
            %        + linearCost(CostMatrix,(bkp1-Xk)));

           fprintf('alpha : %f \n', alpha);
            %fprintf('alpha opt %f : \n',  alpha_opt);

            
            %fprintf('t0 %f : \n', t0);

            %fprintf('beta %f : \n',  beta__opt);
           fprintf('beta opt %f : \n', beta_);

            % 

           fprintf('--------------------------\n');
           pause;
        end
        if ((beta_ < 0.00001) || (t0 >= 1) || ( k +1 == maxIter))
            %fprintf('No line search\n');
            Xkp1 = bkp1;
            S(k+1) = S(k+1); %% Already storedx
            Lterm = Lterm_new;
        else
            
            %fprintf('line search\n');
            Xkp1 = Xk + t0*(bkp1 - Xk); 
            S(k+1) = S(k) - ((alpha^2)/(4*beta_)); %%tester S(k+1)
            Lterm = linearCost(CostMatrix, Xkp1);
        end
        k = k+1;
        Xkm1 = Xk;
        Xk=Xkp1;
        

        %%%%%%%% DEBUG %%%%%%%%%%%
        % S'
        % Lterm
        % alpha
        %beta_
        %fprintf('diff k,k-1 : %d\n',  norm(Xk - Xkm1,1));
        %fprintf('S(k) : %f \n', S(k));
        
        % true_sk = valueObjFunction(G1, G2, CostMatrix, removeEpsilonAssociation(G1,G2,Xk), ...
        %                            cei, ced, ces);
        
        %fprintf('value Obj function computed from Xk  : %f \n', true_sk);

        % if (true_sk != S(k))
        %    fprintf('Problem : true_sk and S(k) should be equal \n');
        %     true_sk
        %     S(k)
        %     Xk
        % end
        
        %fprintf('Condition sur itération : %d \n', ( k <= maxIter));
        %fprintf('Norme de différences des deux matrices : %d \n', (norm(Xk- Xkm1,1) > 0.0001));
        % pause;
        %%%%%%%% END DEBUG %%%%%%%%%%% 
        % Xk
        % bkp1
        % pause;
    end
    
    running_time=toc(chrono);
    mapping = bkp1;
    cost = S(k);
