%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Produce Results for Draft 
%% Results are converted into tables and figures in Draft_Graphs_Tables



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Histogram Iteration Method (Line 168 in Julia)

% Run Histogram Simulation for Different Grids 
H_grid_size = [250 250 500 1000 5000] ; 
n_H = numel(H_grid_size) ;
H_Gamma_timed = zeros(n_H,1)    ; 
H_Gamma_bytes = zeros(n_H,1)    ;
H_M_timed = zeros(n_H,3)  ; 
H_M_bytes = zeros(n_H,3)  ;

H_a_grid  = zeros(H_grid_size(end),n_H) ; % Asset Grid
H_Gamma_a     = zeros(H_grid_size(end),n_H) ; % Distribution of assets

H_Wealth_Stats = zeros(n_H,6) ; % 5 Percentiles (levels) + Average Assets
H_Wealth_Share = zeros(n_H,5) ; % 5 Top Shares by percentiles
H_Pareto_Coeff = zeros(n_H,1  ) ; % Pareto Coefficient
pct_list = [90 95 99 99.9 99.99] ; % Percentiles to be computed

H_Decile       = zeros(11,2 ,n_H) ;
H_Decile_Tr    = zeros(10,10,n_H) ;
N_Decile_Tr    = 9                ;

H_Cons_Corr    = zeros(n_H,1);
H_A_Corr  = zeros(n_H,1) ;
H_eps_Corr = zeros(n_H,1);
H_zeta_Corr  = zeros(n_H,1);
N_Cons_Corr = 2 ;

for i=1:n_H

    disp("                    ")
    fprintf(' Histogram with (%i ) Grid Points \n', H_grid_size(i))
    disp("                    ")

    % Set up model structure
    M_Hist = Functions_ModelSolution.model2modify(M, 'read_flag',1 );
    M_Hist = Functions_ModelSolution.model2modify(M_Hist, 'n_a_fine',H_grid_size(i) );
    
    initial_tic=tic;
    M_Hist = Functions_ModelSolution.Aiyagari_Equilibrium(M_Hist);
    H_Gamma_timed(i) = toc(initial_tic);

    % Grid and distribution
    H_a_grid(1:H_grid_size(i),i) = M_Hist.a_grid_fine ;
    H_Gamma_a(1:H_grid_size(i),i) = squeeze(sum(M_Hist.Gamma, [3,2])); 

    % Moments
    initial_tic=tic;
    [H_Wealth_Stats(i,:), H_Wealth_Share(i,:), H_Pareto_Coeff(i)] = Functions_DraftResults.H_Moments_Top_Shares(M_Hist,pct_list) ; 
    H_M_timed(i,1) = toc(initial_tic);

    % Decile Transitions
    initial_tic=tic;
    [H_Decile_Tr(:,:,i), H_Decile(:,:,i) ] = Functions_DraftResults.H_Moments_Decile_Transitions(M_Hist,N_Decile_Tr) ; 
    H_M_timed(i,2) = toc(initial_tic);
    
    % Consumption Autocorrelation for first quintile
     initial_tic=tic;
    [ H_Cons_Corr(i), H_A_Corr(i), H_eps_Corr(i), H_zeta_Corr(i) ] = Functions_DraftResults.H_Moments_Auto_Correlation(M_Hist,N_Cons_Corr,1,int32(H_Decile(3,1,i))) ;
    H_M_timed(i,3) = toc(initial_tic);

end % end of loop i=1:n_H

writematrix(H_Gamma_timed,'./Hist_Folder/H_G_timed.csv')
writematrix(H_Gamma_bytes,'./Hist_Folder/H_G_bytes.csv')
writematrix(H_a_grid,'./Hist_Folder/H_a_grid.csv')
writematrix(H_Gamma_a,'./Hist_Folder/H_G_a.csv')
writematrix(H_M_timed,'./Hist_Folder/H_M_timed.csv')
writematrix(H_M_bytes,'./Hist_Folder/H_M_bytes.csv')
writematrix(H_Wealth_Stats,'./Hist_Folder/H_Wealth_Stats.csv')
writematrix(H_Wealth_Share,'./Hist_Folder/H_Wealth_Share.csv')
writematrix(H_Pareto_Coeff,'./Hist_Folder/H_Pareto_Coeff.csv')
writematrix(H_Decile,'./Hist_Folder/H_Decile.csv')
writematrix(H_Decile_Tr,'./Hist_Folder/H_Decile_Tr.csv')
writematrix(H_Cons_Corr,'./Hist_Folder/H_Cons_Corr.csv')
writematrix(H_A_Corr,'./Hist_Folder/H_A_Corr.csv')
writematrix(H_eps_Corr,'./Hist_Folder/H_eps_Corr.csv')
writematrix(H_zeta_Corr,'./Hist_Folder/H_z_Corr.csv')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Monte Carlo Simulation 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run Simulation for Different Panel Size (Line 300 in Julia)

    % Set up model structures 
    M_Simul = Functions_ModelSolution.ModelStructure(par, ...
    120000, 4.5, 4.5, 250, 500, 7,11,1,false);
    M_Simul = Functions_ModelSolution.model2modify(M_Simul, 'read_flag',1 );
    % M_Panel = Model_Panel(N_Panel=1000000)   ; 

    % Set up discrete observations 
    S_sample = [50000 250000 500000 1000000 10000000]    ; 
    N_S      = numel(S_sample)         ;
    pct_list = [90;95;99;99.9;99.99]    ;  

    S_M_timed      = zeros(N_S,4)       ; % 1-> Simulation 2->Top Shares 3->Decile Transition 4->Auto-correlation
    S_M_bytes      = zeros(N_S,4)       ; % 1-> Simulation 2->Top Shares 3->Decile Transition 4->Auto-correlation
    
    S_Wealth_Sample= zeros(N_S,S_sample(N_S)) ;
    S_Wealth_Stats = zeros(N_S,6)       ; 
    S_Wealth_Share = zeros(N_S,5)       ; 
    S_Pareto_Coeff = zeros(N_S  )       ;  

    S_Decile       = zeros(11,N_S)      ;
    S_Decile_Tr    = zeros(10,10,N_S)   ;

    S_Cons_Corr    = zeros(N_S)         ;
    S_A_Corr       = zeros(N_S)         ;
    S_eps_Corr       = zeros(N_S)         ;
    S_zeta_Corr       = zeros(N_S)         ;
    
    % Solve model 
    initial_tic=tic;
    M_Simul = Functions_ModelSolution.Aiyagari_Equilibrium(M_Simul);    
    S_Gamma_timed=toc(initial_tic);



    % Moments
    for i=1:N_S

        disp(" ")
        fprintf(' Simulation with %i agents \n', S_sample(i))
        disp(" ")

    % Simulate Panel
        M_Panel = Functions_MonteCarlo.PanelStructure(500000, 10, 1000, 1000, 1e-2, 3489398);
        M_Panel = Functions_MonteCarlo.panelmodel2modify(M_Panel,'N_Panel',S_sample(i));
        initial_tic=tic;
        M_Panel = Functions_MonteCarlo.Simulate_Panel(M_Simul,M_Panel, 0);
        S_M_timed(i,1)=toc(initial_tic);

        %S_M_bytes(i,1) = 0;
        

    % 1-2) Top Wealth Shares and Pareto Coefficient 
        initial_tic=tic;
        

        % Select sample 
        a_sample  = M_Panel.a_mat(:,end) ;% M_Panel.a_mat[1:S_sample[i],end] ;
        S_Wealth_Sample(i,1:S_sample(i)) = a_sample  ;

        % Average wealth 
        S_Wealth_Stats(i,end)     = mean( a_sample )   ;

        % Percentiles 
        S_Wealth_Stats(i,1:end-1) = prctile( a_sample , pct_list ) ;

        % Top Shares 
        for p=1:5
            S_Wealth_Share(i,p)= 100*sum( a_sample(a_sample>=S_Wealth_Stats(i,p))) /(S_sample(i)*S_Wealth_Stats(i,end));
        end
        

        % Pareto Coefficient 
        Pareto_sample = sort( a_sample(a_sample>=1000 ) ) ; % Select and sort observations above $1M  
        Pareto_CCDF   =  transpose((numel(Pareto_sample):-1:1)/numel(Pareto_sample)); % Counter CDF = 1- CDF
        S_Pareto_Coeff(i) = dot(log(Pareto_sample(1:end-1)/1000),log(Pareto_sample(1:end-1)/1000)) \ dot(log(Pareto_sample(1:end-1)/1000),log(Pareto_CCDF(1:end-1))) ;

        % Time it 
        S_M_timed(i,2) = toc(initial_tic);
        %S_M_bytes(i,2)




    % 3) Decile Transitions 
        initial_tic=tic;
        
        a_aux_0       = M_Panel.a_mat(:,1)   ;
        a_aux_T       = M_Panel.a_mat(:,end) ;
        S_Decile(:,i) = prctile( a_aux_T , 0:10:100 ) ; % Deciles based on end of sample (hoping for stationarity)
        for p=1:10
            ind_d = find( (S_Decile(p,i)<=a_aux_0 ) & (a_aux_0 <= S_Decile(p+1,i)) ); % Find "i" in each decile in t=1
            a_T   = a_aux_T(ind_d) ; % Follow them to t=10
            for pT=1:10
                inner_ind = (S_Decile(pT,i)<=a_T) & (a_T <= S_Decile(pT+1,i)) ;
                S_Decile_Tr(p,pT,i) = 100*sum(inner_ind)/numel(ind_d); % Get transition rates 
            end

        end 

        % Time it 
        S_M_timed(i,3) = toc(initial_tic);
        %S_M_bytes(i,3)

    % 4) Consumption Autocorrelation for first quintile
        initial_tic=tic;
        % Fix current sample and future sample 
        a_aux_0 = M_Panel.a_mat(:, 8 ) ; 
        c_aux_0 = M_Panel.c_mat(:, 8 ) ; 
        eps_aux_0 = M_Simul.eps_grid(M_Panel.eps_mat(:, 8 )) ; 
        zeta_aux_0 = M_Simul.zeta_grid(M_Panel.zeta_mat(:, 8 ))   ;
        a_aux_T = M_Panel.a_mat(:,end) ; 
        c_aux_T = M_Panel.c_mat(:,end) ; 
        eps_aux_T = M_Simul.eps_grid(M_Panel.eps_mat(:,end)) ; 
        zeta_aux_T = M_Simul.zeta_grid(M_Panel.zeta_mat(:,end)) ;

        % Find index of first quintile 
        ind_q  = find( (0 <= a_aux_0 )  & (a_aux_0<=S_Decile(3,i) )); % Index of first quintile of assets in T-1
        % Compute moments 
        inner = corrcoef(c_aux_0(ind_q) , c_aux_T(ind_q))  ;
        S_Cons_Corr(i)= inner(1,2) ;
        inner = corrcoef( a_aux_0(ind_q), a_aux_T(ind_q) );
        S_A_Corr(i)   = inner(1,2) ;
        inner = corrcoef( log(eps_aux_0(ind_q)), log(eps_aux_T(ind_q))) ;
        S_eps_Corr(i)   = inner(1,2) ;
        inner = corrcoef( log(zeta_aux_0(ind_q)), log(zeta_aux_T(ind_q))) ;
        S_zeta_Corr(i)   = inner(1,2) ; 
        clear inner 

        % Deconstruct measure: mean, var, cov 
        av_c_0   = mean(     c_aux_0(ind_q) ) ; 
        av_c_T   = mean(     c_aux_T(ind_q) ) ; 
        cov_c_0T = cov(c_aux_0(ind_q), c_aux_T(ind_q)) ;  
        av_a_0   = mean(     a_aux_0(ind_q) ) ; 
        av_a_T   = mean(     a_aux_T(ind_q) ) ; 
        cov_a_0T = cov(a_aux_0(ind_q), a_aux_T(ind_q)) ;  
        av_eps_0   = mean(log(eps_aux_0(ind_q))) ; 
        av_eps_T   = mean(log(eps_aux_T(ind_q))) ; 
        cov_eps_0T = cov(log(eps_aux_0(ind_q)), log(eps_aux_T(ind_q)) )  ;  
        av_zeta_0   = mean(log(zeta_aux_0(ind_q))) ; 
        av_zeta_T   = mean(log(zeta_aux_T(ind_q))) ; 
        cov_zeta_0T = cov((log(zeta_aux_0(ind_q))), (log(zeta_aux_T(ind_q)))) ;  
        disp("Moments for autocorrelation of consumption:")
        fprintf('cor_c= %3.2f - av_c_0 = %5.3f - av_c_T= %5.3f - sd_c_0= %3.2f  - sd_c_T=  %4.2f   - cov_0N= %5.3f   \n', ...
            round(100*S_Cons_Corr(i),2 ), round(av_c_0,3), round(av_c_T,3), round(sqrt(cov_c_0T(1,1)), 2), ...
            round(sqrt(cov_c_0T(2,2)),2), round(cov_c_0T(1,2),3) )
          
        fprintf('cor_a= %3.2f - av_a_0 = %5.3f - av_a_T= %5.3f - sd_a_0= %3.2f  - sd_a_T=  %4.2f   - cov_0N= %5.3f   \n', ...
            round(100*S_A_Corr(i),2 ), round(av_a_0,3), round(av_a_T,3), round(sqrt(cov_a_0T(1,1)), 2), ...
            round(sqrt(cov_a_0T(2,2)),2), round(cov_a_0T(1,2),3) )
          

        fprintf('cor_eps= %3.2f - av_eps_0 = %5.3f - av_eps_T= %5.3f - sd_eps_0= %3.2f  - sd_eps_T=  %4.2f   - cov_0N= %5.3f   \n', ...
            round(100*S_eps_Corr(i),2 ), round(av_eps_0,3), round(av_eps_T,3), round(sqrt(cov_eps_0T(1,1)), 2), ...
            round(sqrt(cov_eps_0T(2,2)),2), round(cov_eps_0T(1,2),3) )
          

        fprintf('cor_zeta= %3.2f - av_zeta_0 = %5.3f - av_zeta_T= %5.3f - sd_zeta_0= %3.2f  - sd_zeta_T=  %4.2f   - cov_0N= %5.3f   \n', ...
            round(100*S_zeta_Corr(i),2 ), round(av_zeta_0,3), round(av_zeta_T,3), round(sqrt(cov_zeta_0T(1,1)), 2), ...
            round(sqrt(cov_zeta_0T(2,2)),2), round(cov_zeta_0T(1,2),3) )
          
        S_M_timed(i,4)=  toc(initial_tic);
    end % End of loop  i=1:N_S



    writematrix(S_M_timed,'./MC_Folder/S_M_timed.csv')
    %writematrix(S_M_bytes,'./MC_Folder/S_M_bytes.csv')
    writematrix(S_Gamma_timed ,'./MC_Folder/S_G_timed.csv')
    %writematrix(S_Gamma_bytes,'./MC_Folder/S_G_bytes.csv')





    S_Wealth_1 = S_Wealth_Sample(1:N_S-1,1:S_sample(N_S-1)) ; 
    S_Wealth_2 = S_Wealth_Sample(N_S,:);

    
    writematrix(S_Wealth_1,'./MC_Folder/S_Wealth_Sample.csv')
    
    writematrix(round(S_Wealth_2,2),'./MC_Folder/S_Wealth_Sample_Large.csv')

    writematrix(S_Wealth_Stats,'./MC_Folder/S_Wealth_Stats.csv')

    writematrix(S_Wealth_Share,'./MC_Folder/S_Wealth_Share.csv')

    writematrix(S_Pareto_Coeff,'./MC_Folder/S_Pareto_Coeff.csv')

    writematrix(S_Decile,'./MC_Folder/S_Decile.csv')

    writematrix(S_Decile_Tr,'./MC_Folder/S_Decile_Tr.csv')
    
    writematrix(S_Cons_Corr,'./MC_Folder/S_Cons_Corr.csv')
    
    writematrix(S_A_Corr,'./MC_Folder/S_A_Corr.csv')
        
    writematrix(S_eps_Corr,'./MC_Folder/S_eps_Corr.csv')

    writematrix(S_zeta_Corr,'./MC_Folder/S_z_Corr.csv')

   




