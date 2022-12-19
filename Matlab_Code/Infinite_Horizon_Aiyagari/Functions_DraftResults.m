
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions for Monte-Carlo simulation 


    classdef Functions_DraftResults
        methods(Static)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            function [Wealth_Stats, Wealth_Share, P_coeff] = H_Moments_Top_Shares(Model, pct_list)
                % Local variable
                M_in = Model;
                % ---------- Output ----------
                % Wealth_Stats, Wealth_Share, P_coeff    


                % Allocate results
                 Wealth_Stats = zeros(numel(pct_list)+1,1) ;
                 Wealth_Share = zeros(numel(pct_list),1) ;

                 % Compute distribution of assets
                 Gamma_a = squeeze(sum(M_in.Gamma, [3,2]));
                 CDF_a = cumsum(Gamma_a);

                 % Compute moments from distribution
                 Wealth_Stats(end) = dot(M_in.a_grid_fine, Gamma_a);
                 for j=1:numel(pct_list)
                     ind_top = (CDF_a>=pct_list(j)/100);
                     WS = transpose(M_in.a_grid_fine(ind_top));
                     Wealth_Stats(j) = WS(1);
                     Wealth_Share(j) = 100*dot(M_in.a_grid_fine(ind_top),Gamma_a(ind_top))/Wealth_Stats(end) ;
                     clear WS
                 end

                 % Pareto Coefficient
                 ind = transpose(find(M_in.a_grid_fine>=1000));
                 grid_1M = M_in.a_grid_fine(ind) ;
                 Gamma_a_1M = Gamma_a(ind)/sum(Gamma_a(ind));
                 Gamma_a_1M = Gamma_a_1M/sum(Gamma_a_1M);
                 CCDF_1M = 1 - cumsum(Gamma_a_1M);
                 P_coeff = dot(log(grid_1M(1:end-1)/1000),log(grid_1M(1:end-1)/1000)) \ dot(log(grid_1M(1:end-1)/1000),log(CCDF_1M(1:end-1)) ) ;
                 
                 % Return Moments
                 return
            
            end % End of function H_Moments_Top_Shares



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Tr_deciles_a, deciles_mat] = H_Moments_Decile_Transitions(Model, N)
    
    % Local Variable
    M_in=Model;
    
    % ---------- Output ----------
    % Output 1 Tr_deciles_a
    % Output 2 deciles_mat

    % Allocate results
    Tr_deciles_a = zeros(10,10) ;
    
    % Compute distribution of assets 
    Gamma_a = squeeze(sum(M_in.Gamma, [3,2]));
    CDF_a = cumsum(Gamma_a);

    % Find deciles
    deciles_a = randi(10,[11,1]);
    deciles_l = randi(10,[11,1]);
    deciles_a(1) = 0;


    for i=2:10
        ind = find((100*CDF_a)>=(10*(i-1)));
        collected = transpose(1:1:M_in.n_a_fine);
        deciles_a(i) = collected(ind(1));
        deciles_l(i) = M_in.a_grid_fine(deciles_a(i));
        clear ind 
    end
    deciles_a(11)=M_in.n_a_fine;
    deciles_l(11)=M_in.a_grid_fine(end);
    deciles_mat = horzcat(deciles_a, deciles_l);
    
   
    %  For each decile compute transitions   
    for i=1:10
        % Conditional distribution on current decile 
            Gamma_0 = zeros(M_in.n_a_fine, M_in.n_eps, M_in.n_zeta);
            Gamma_0(deciles_a(i)+1:deciles_a(i+1), :, :) = M_in.Gamma(deciles_a(i)+1:deciles_a(i+1), :, :);
            Gamma_0 = Gamma_0./sum(sum(sum(Gamma_0)));
            % Iterate distribution
            Gamma_N = Functions_ModelSolution.Histogram_Iteration(M_in, N, Gamma_0);
            % Obtain transitions (sum over new distribution within deciles)
            inner_column = zeros(10,1);
            for j=1:10
                inner_column(j) = 100*sum(sum(sum(Gamma_N(deciles_a(j)+1: deciles_a(j+1),:,:))));
            end
            Tr_deciles_a(i,:) = inner_column;
            clear inner_column;
    end
    

end % End of function H_Moments_Decile_Transitions



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cor_c_q, cor_a_q, cor_eps_q, cor_zeta_q] = H_Moments_Auto_Correlation(Model,N,a_min,a_max)
    % Local Variable
    M_in = Model;

    % Get average consumption and standard deviation
    % Initial distribution over asset range
    Gamma_q = zeros(M_in.n_a_fine, M_in.n_eps, M_in.n_zeta);
    Gamma_q(a_min:a_max, :,:) = M_in.Gamma(a_min:a_max,:,:);
    Gamma_q = Gamma_q/sum(sum(sum(Gamma_q)));

    % Consumption
    av_c_q = sum( sum( sum( M_in.G_c_fine(a_min:a_max,:,:).*Gamma_q(a_min:a_max,:,:) ) ));
    sd_c_q = sqrt( sum(sum(sum( ((M_in.G_c_fine(a_min:a_max,:,:) - av_c_q).^2).*Gamma_q(a_min:a_max,:,:) ) ))) ;
    % Assets
    av_a_q = sum(sum(sum( M_in.a_mat_fine(a_min:a_max,:,:).*Gamma_q(a_min:a_max,:,:) ))) ;
    sd_a_q = sqrt( sum(sum(sum( ((M_in.a_mat_fine(a_min:a_max,:,:)-av_a_q).^2).*Gamma_q(a_min:a_max,:,:) ) ))) ;
    % Epsilon
    av_eps_q = sum(sum(sum( log(M_in.eps_mat_fine(a_min:a_max,:,:)).*Gamma_q(a_min:a_max,:,:) ))) ;
    sd_eps_q = sqrt(sum( sum(sum( ((log(M_in.eps_mat_fine(a_min:a_max,:,:)) -av_eps_q).^2).*Gamma_q(a_min:a_max,:,:) ) ))) ;
    % Zeta
    av_zeta_q = sum(sum(sum(log(M_in.zeta_mat_fine(a_min:a_max,:,:)).*Gamma_q(a_min:a_max,:,:) ))) ;
    sd_zeta_q = sqrt( sum(sum(sum( ((log(M_in.zeta_mat_fine(a_min:a_max,:,:)) -av_zeta_q).^2).*Gamma_q(a_min:a_max,:,:) ) ))) ;

    % Future average and standard deviation conditional on first quintile
    % Iterate distribution 

    Gamma_qN = Functions_ModelSolution.Histogram_Iteration(M_in,N, Gamma_q);
    % Consumption
     av_c_q_N = sum(sum(sum( M_in.G_c_fine.*Gamma_qN  )));
     sd_c_q_N = sqrt( sum(sum(sum( ((M_in.G_c_fine-av_c_q_N).^2).*Gamma_qN  ))) )   ; 
     % Assets 
     av_a_q_N = sum(sum(sum( M_in.a_mat_fine.*Gamma_qN  )));
     sd_a_q_N = sqrt( sum(sum(sum( ((M_in.a_mat_fine-av_a_q_N).^2).*Gamma_qN  ) ) )); 
     % epsilon 
     av_eps_q_N = sum(sum(sum( log(M_in.eps_mat_fine).*Gamma_qN  )));
     sd_eps_q_N = sqrt(sum(sum( sum( ((log(M_in.eps_mat_fine)-av_eps_q_N).^2).*Gamma_qN  ) ))); 
     % zeta
     av_zeta_q_N = sum(sum(sum( log(M_in.zeta_mat_fine).*Gamma_qN  )));
     sd_zeta_q_N = sqrt( sum(sum(sum( ((log(M_in.zeta_mat_fine)-av_zeta_q_N).^2).*Gamma_qN  ) ))); 

     % Compute integrand of correlation
    cov_c = zeros(M_in.n_a_fine,M_in.n_eps,M_in.n_zeta) ; 
    cov_a = zeros(M_in.n_a_fine,M_in.n_eps,M_in.n_zeta) ; 
    cov_eps = zeros(M_in.n_a_fine,M_in.n_eps,M_in.n_zeta) ; 
    cov_zeta = zeros(M_in.n_a_fine,M_in.n_eps,M_in.n_zeta) ; 

    for i_zeta =1:M_in.n_zeta % Current zeta
        for i_eps = 1:M_in.n_eps % Current epsilon
        for i_a =a_min:a_max % Current a
            % For each state in state space (a,epsilon,zeta) get conditional distribution 
            Gamma_0 = zeros(M_in.n_a_fine,M_in.n_eps,M_in.n_zeta) ; 
            Gamma_0(i_a,i_eps,i_zeta) = 1 ; 
            % Iterate distribution
            Gamma_N = Functions_ModelSolution.Histogram_Iteration(M_in,N, Gamma_0 ) ;
            % Get portion of integrand (C_0 - av_c)*(C_n - av_c)*Gamma_n
            cov_c(i_a,i_eps,i_zeta) = sum(sum(sum( (      M_in.G_c_fine(i_a,i_eps,i_zeta)   -av_c_q )*(      M_in.G_c_fine   -av_c_q_N ).*Gamma_N ))) ;
            cov_a(i_a,i_eps,i_zeta) =sum(sum( sum( (      M_in.a_mat_fine(i_a,i_eps,i_zeta) -av_a_q )*(      M_in.a_mat_fine -av_a_q_N ).*Gamma_N ))) ;
            cov_eps(i_a,i_eps,i_zeta) = sum(sum(sum( ( log(M_in.eps_mat_fine(i_a,i_eps,i_zeta))-av_eps_q )*( log(M_in.eps_mat_fine)-av_eps_q_N ).*Gamma_N ))) ;
            cov_zeta(i_a,i_eps,i_zeta) = sum(sum(sum( ( log(M_in.zeta_mat_fine(i_a,i_eps,i_zeta))-av_zeta_q )*( log(M_in.zeta_mat_fine)-av_zeta_q_N ).*Gamma_N ))) ;
        end % End of loop a
        end % End of loop eps
    end % End of loop zeta


    % Integrate with respect to initial distribution
    cov_c_q  = sum(sum(sum( cov_c.*Gamma_q ))) ;   
    cor_c_q  = 100*cov_c_q/(sd_c_q*sd_c_q_N) ; 
    cov_a_q  =  sum(sum(sum( cov_a.*Gamma_q ) )) ;    
    cor_a_q  = 100*cov_a_q/(sd_a_q*sd_a_q_N) ; 
    cov_eps_q  =  sum(sum(sum( cov_eps.*Gamma_q )))  ;    
    cor_eps_q  = 100*cov_eps_q/(sd_eps_q*sd_eps_q_N) ; 
    cov_zeta_q  =  sum(sum(sum( cov_zeta.*Gamma_q )))  ;  
    cor_zeta_q  = 100*cov_zeta_q/(sd_zeta_q*sd_zeta_q_N) ; 
    
    % Print other moments
    fprintf('cor_c = %4.2f - cor_a = %4.2f - cor_eps = %4.2f - cor_zeta = %4.2f - a in [%8.2f %8.2f] \n ', ...
        round(cor_c_q,2), round(cor_a_q,2), round(cor_eps_q,2), round(cor_zeta_q,2), round(M_in.a_mat_fine(a_min),2), round(M_in.a_mat_fine(a_max),2) )

    
    

end % End of function H_Moments_Auto_Correlation


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% End of file
        end % End of methods(Static)
    end %End of Functions_DraftResults
