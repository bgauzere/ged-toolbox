function [mapping, running_time,conv,S_zeta] = ipfp_zeta(G1, G2, ...
                                                      costs, maxIter, Mapping_init,zeta, debug)
    %epsilon = - sign(zeta);
    nu = abs(zeta);
    chrono=tic;
    if nargin < 7
        debug = 0;
    end
    R_zeta= zeros(maxIter,1);
    S = zeros(maxIter,1);
    S_zeta = zeros(maxIter,1);
    
    tic;
    % n=size(G1,1);
    % m=size(G2,1);
    CostMatrix = NodeCostMatrix(G1, G2, costs);
    Xk = Mapping_init; 
    [phi_i, i] = find(Xk');
    % Xk
    XkD = quadraticTerm(G1, G2, costs.cei, costs.ced, costs.ces, Xk, ...
                        [i, phi_i]);
    Lterm = linearMul(CostMatrix, Xk);
    S(1) = 0.5*linearMul(XkD,Xk) + Lterm; % Inutile
    S_zeta(1) = (1-nu)*S(1) + zeta*linearMul(Xk,Xk); % Premier terme inutile
    % S_zeta(1) = zeta*linearMul(Xk,Xk);
    k=1;
    Xkm1 = 2*Xk; % For the first loop
    %% !! Premiere iter, 1 + epsilon*nu = 0 => mapping init est
    %% inutile, possibilité de gagner du temps ?
    while ((k < maxIter)  && (norm(Xk - Xkm1,1) > 0.0001));
        [phi_i, i] = find(Xk');
        XtD = quadraticTerm(G1, G2, costs.cei, costs.ced, costs.ces, Xk, ...
                            [i, phi_i]);
        linear_sub_problem = (1-nu) * (XtD + CostMatrix) + (2 * zeta*Xk); 
        % if (debug)
        %     linear_sub_problem
        % end
        lsap_sub = linear_sub_problem;
        %%Translation to a positive cost matrix
        min_value = min(min(linear_sub_problem)); 
        lsap_sub = lsap_sub + abs(min_value);
        lsap_sub(CostMatrix > 1000) = -1;

        [phi_sub_problem,~,~] = hungarianLSAP(lsap_sub);
        nplusm = length(phi_sub_problem);
        bkp1=zeros(nplusm,nplusm);
        %%bkp1 is an assignment matrix
        %Computation of assignement matrix from phi
        bkp1(sub2ind([nplusm,nplusm],int32((1:nplusm)'),phi_sub_problem+1)) = 1;%Computation of assignement matrix from phi
        % R(k+1) = linearMul((XtD + CostMatrix),bkp1);
   	R_zeta(k+1) = linearMul(linear_sub_problem,bkp1);
        
        Lterm_new  = linearMul(CostMatrix, bkp1);
        [phi_i, i] = find(bkp1');
        Xkp1tD = quadraticTerm(G1, G2, costs.cei, costs.ced, costs.ces, bkp1,[i, phi_i]);
        quadratic_part = linearMul( Xkp1tD,bkp1);
        S(k+1) =  0.5*quadratic_part + Lterm_new;
        % S_bkp1 = S(k+1);
        S_zeta(k+1) = (1-nu)*S(k+1) + zeta*linearMul(bkp1,bkp1);
        % S_zeta_bkp1 = S_zeta(k+1);
        % current_ged = quadratic_part/2 + Lterm_new;
        
        
        alpha_zeta = R_zeta(k+1) - 2*S_zeta(k) + (1-nu)*Lterm;
        beta_zeta = S_zeta(k+1) + S_zeta(k) - R_zeta(k+1) - (1-nu)* Lterm;
        t0 = 0;
        if (beta_zeta > 0.00001)
            %% Compute t0 here to avoid division by 0 
            t0 = -alpha_zeta / (2* beta_zeta); 
        end 
        % S_bkp1 = S(k+1); 
        if ((beta_zeta < 0.00001) || (t0 >= 1))
            Xkp1 = bkp1;
            % S_zeta(k+1) = S_zeta(k+1); Useless
            % if (debug)
            %     LtermXk = Lterm;
            % end
            Lterm = Lterm_new;
        else %% line search
            Xkp1 = Xk + t0*(bkp1 - Xk); 
            % if (debug)
            %     [phi_i, i] = find(Xkp1');
            %     tmp = quadraticTerm(G1, G2, costs.cei, costs.ced, costs.ces, Xkp1, [i, phi_i]);
            %     tmp2 = linearMul(tmp,Xkp1);
            %     S(k+1) = 0.5*tmp2 + linearMul(CostMatrix, Xkp1);
            % end
            S_zeta(k+1) = S_zeta(k) - ((alpha_zeta^2)/(4*beta_zeta));
            % LtermXk = Lterm;
            Lterm = linearMul(CostMatrix, Xkp1);
        end
        % if (debug)
        %     fprintf('zeta = %f \n',zeta);
        %     fprintf('1 - nu = %f \n',1 - nu);
        %     fprintf('R_zeta(xk)  = %f \n',linearMul(linear_sub_problem,Xk));
        %     fprintf('R_zeta(k+1) = %f \n',R_zeta(k+1));

        %     % fprintf('S(k) = %f \n',S(k));
        %     % fprintf('S(k+1) = %f  \n', S(k+1));

        %     fprintf('S_zeta(k)    = %f \n',S_zeta(k));
        %     fprintf('S_zeta(k+1)  = %f \n', S_zeta(k+1));
            
        %     % fprintf('Lterm = %f \n',Lterm);
	%     % fprintf('Quadratic part = %f \n',quadratic_part);
        %     % fprintf('my current ged = %f \n', sum(sum( Xkp1tD.*bkp1))/2 + linearCost(CostMatrix, bkp1));
            
        %     alpha = R(k+1) - 2*S(k) + LtermXk;
        %     [phi_i, i] = find(Xk');
        %     tmp = quadraticTerm(G1, G2, costs.cei, costs.ced, costs.ces, Xk, [i, phi_i]);
        %     tmp2 = linearMul( tmp,Xk);
        %     alpha = linearMul(XtD + CostMatrix, bkp1) - 2* (0.5*tmp2 + linearMul(CostMatrix,Xk)) + linearMul(CostMatrix, Xk);
            
        
        %     [phi_i, i] = find(bkp1');
        %     tmp_bkp1 = quadraticTerm(G1, G2, costs.cei, costs.ced, costs.ces, bkp1, ...
        %                         [i, phi_i]);
        %     tmp2_bkp1 = linearMul( tmp_bkp1,bkp1);
            
        %     fprintf('R_zeta(k+1) - R_zeta(xk) = %f \n',R_zeta(k+1) ...
        %             - linearMul(linear_sub_problem,Xk));
        %     fprintf('R_zeta(bk+1 -xk)         = %f \n', ...
        %             linearMul(linear_sub_problem,bkp1 - Xk));
        %     fprintf('alpha_zeta               = %f \n', ...
        %             alpha_zeta);
        %     alpha_zeta_bis= (1-nu)*alpha + 2*zeta*linearMul(Xk,(bkp1 - Xk));
        %     fprintf('alpha_zeta (bis)         = %f \n', ...
        %             alpha_zeta_bis);
            
        %     alpha_zeta_ter = linearMul(linear_sub_problem, bkp1) - 2 * ...
        %         ((1-nu)*(0.5*tmp2 + linearMul(CostMatrix,Xk)) + zeta ...
        %          * linearMul(Xk,Xk)) + (1-nu)*linearMul(CostMatrix, ...
        %                                                 Xk);
        %     fprintf('alpha_zeta (ter)         = %f \n', alpha_zeta_ter);
        %     % fprintf('alpha_zeta (bis) =  %f + %f \n', (1-nu)*alpha ,2*zeta*sum(sum(Xk.*(bkp1 - Xk))));
        %     % fprintf('alpha_zeta = %f -%f + %f \n', R_zeta(k+1), 2*S_zeta(k),(1-nu)*Lterm);
            
        %     fprintf('beta_zeta       = %f  \n', beta_zeta);
        %     beta_ = (0.5*tmp2_bkp1  - linearMul(tmp,bkp1) + 0.5*tmp2);
        %     beta_zeta_bis = (1 - nu)*beta_ + zeta*(linearMul((bkp1- Xk),(bkp1 - Xk)));
        %     fprintf('beta_zeta (bis) = %f  \n',beta_zeta_bis);
            
        %     fprintf('t0       = %f \n', t0);
        %     t0_bis =  - alpha_zeta_bis / (2* beta_zeta_bis);
        %     fprintf('t0 (bis) = %f \n', t0_bis);
            

        %     Xkp1_bis = Xk + t0_bis*(bkp1 - Xk); 
        %     [phi_i, i] = find(Xkp1_bis');
        %     Xkp1D = quadraticTerm(G1, G2, costs.cei, costs.ced, costs.ces, Xkp1_bis, [i, phi_i]);
        %     quadratic_part = linearMul(Xkp1D,Xkp1_bis);
        %     Lterm_new_bis = linearMul(CostMatrix,Xkp1_bis);
        %     S_bis = 0.5*quadratic_part + Lterm_new_bis;
        %     S_zeta_bis = (1-nu)*S_bis + zeta*linearMul(Xkp1_bis,Xkp1_bis);
        %     fprintf('S_zeta(k+1)          = %f  \n', S_zeta(k+1));
        %     fprintf('S_zeta(k+1) (bis)    = %f  \n', S_zeta_bis);
        %     S_zeta_bis = S_zeta(k) - ((alpha_zeta_bis^2)/(4*beta_zeta_bis));
        %     fprintf('S_zeta(k+1) (t0 bis) = %f  \n', S_zeta_bis);
        %     fprintf('S_zeta(k+1) (ter)    = %f  \n', S_zeta(k) + ...
        %             t0*alpha_zeta + (t0^2)*beta_zeta);
            
        %     fprintf('S_zeta(bkp1)    = %f  \n', S_zeta_bkp1);
        %     S_bkp1 = S(k+1);
        %     fprintf('S(bkp1)    = %f  \n', S_bkp1);
        %     fprintf('S(k)    = %f  \n', S(k));
        %     fprintf('S(k+1)    = %f  \n', S(k+1));

                        
        %     % fprintf('S_zeta(k+1) = %f - %f  \n', S_zeta(k) , ((alpha_zeta^2)/(4*beta_zeta)));
        %     % fprintf('S_zeta(k+1) (bis) =  %f*%f + %f  \n', (1-nu), S_bis, zeta*sum(sum(Xkp1.*Xkp1)));
        %     % fprintf('S_zeta(k+1) (bis) =  %f + %f  \n', (1-nu)* S_bis , zeta*sum(sum(Xkp1.*Xkp1)));
        %     fprintf('norm conv : %f \n',norm(Xk - Xkm1,1));
        %     fprintf('--------------------------\n');
        %     pause;
        % end
        k = k+1;
        Xkm1 = Xk;
        Xk=Xkp1;
    end   
    running_time=toc(chrono);
    mapping = Xkp1; 
    % cost = S_zeta(k);
    conv = norm(Xk - Xkm1,1);
