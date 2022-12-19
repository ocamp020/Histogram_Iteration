%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Draft Results File
%%% Produce Results for Draft 
%%% Results are converted into tables and figures in Draft_Graphs_Tables

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Run HIstogram Simulation for Different Grids

H_grid_size =  [250 250 500 1000];
n_H = numel(H_grid_size);
H_Gamma_timed = zeros(n_H,1);
H_Gamma_bytes = zeros(n_H,1);
H_M_timed = zeros(n_H,3);
H_M_bytes = zeros(n_H,3);

pct_list =   [90 95 99 99.9 99.99]' ; 
age_0_Wealth_Profile   =  26       ; 
age_0_Wealth_Corr      =  16       ; 
age_T_Wealth_Corr_low  =  21       ; 
age_T_Wealth_Corr_high =  31       ;


H_Wealth_Profile_NB = zeros(M_Aiyagari.p.Max_Age,6,n_H) ; 
H_Wealth_Profile_45 = zeros(M_Aiyagari.p.Max_Age-age_0_Wealth_Profile+1,6,n_H) ; 
H_Wealth_Corr       = zeros(n_H,2)   ; 


for i=1:n_H
    disp("                            ")
    fprintf(' Histogram with %i Grid Points \n', H_grid_size(i))
    disp("                            ")
    
    % Set up model structure
    M_Hist = M_Aiyagari;
    M_Hist.n_a_fine = H_grid_size(i);
    % Solve for stationary distribution and save time and allocation - Adjust grid size



        %% Moments 
        
        % 1) Wealth Profiles 

        % 2) Wealth Autocorrelation 35-45

        % 3) Wealth Autocorrelation 35-55

end

% Writing matrices
writematrix(H_Gamma_timed, 'Hist_Folder/H_G_timed.csv')
writematrix(H_Gamma_bytes, 'Hist_Folder/H_G_bytes.csv')
writematrix(H_M_timed, 'Hist_Folder/H_M_timed.csv')
writematrix(H_M_bytes, 'Hist_Folder/H_M_bytes.csv')
writematrix(H_Wealth_Profile_NB, 'Hist_Folder/H_Wealth_Profile_NB.csv')
writematrix(H_Wealth_Profile_45, 'Hist_Folder/H_Wealth_Profile_45.csv')
writematrix(H_Wealth_Corr, 'Hist_Folder/H_Wealth_Corr.csv')






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Moments from Monte Carlo Simulation 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run Simulation for Different Panel Size 

% Set up age limits 
age_0_Wealth_Profile   =  26           ; 
age_0_Wealth_Corr      =  16           ; 
age_T_Wealth_Corr_low  =  21           ; 
age_T_Wealth_Corr_high =  31           ; 

% Set up model structures
M_Simul = Model;
M_Panel = Functions_MonteCarlo.Model_Panel();
M_C_45 = struct('N_Panel', 500000, 'T_Panel', M_Simul.p.Max_Age-(age_0_Wealth_Profile-1), 'rng_seed', 297835398);
M_C_45.a_mat    = zeros([M_C_45.N_Panel, M_C_45.T_Panel]);
M_C_45.c_mat    =  zeros([M_C_45.N_Panel, M_C_45.T_Panel]);
M_C_45.eps_mat  = zeros([M_C_45.N_Panel, M_C_45.T_Panel]);
M_C_45.h_mat    =  zeros([M_C_45.N_Panel, M_C_45.T_Panel]);
M_C_45.t_mat    = zeros([M_C_45.N_Panel,1]);


M_C_35 = struct('N_Panel', 500000, 'T_Panel', M_Simul.p.Max_Age-(age_0_Wealth_Corr-1), 'rng_seed', 297835398);
M_C_35.a_mat    = zeros([M_C_35.N_Panel, M_C_35.T_Panel]);
M_C_35.c_mat    =  zeros([M_C_35.N_Panel, M_C_35.T_Panel]);
M_C_35.eps_mat  = zeros([M_C_35.N_Panel, M_C_35.T_Panel]);
M_C_35.h_mat    =  zeros([M_C_35.N_Panel, M_C_35.T_Panel]);
M_C_35.t_mat    = zeros([M_C_35.N_Panel,1]);
            

% Set up discrete observations 
S_sample = [50000 100000 250000 500000] ;   
N_S      = numel(S_sample)              ;
pct_list = [90;95;99;99.9;99.99]         ; 
med_eps    = convert(Int64,round(M.n_eps/2)) ;

S_M_timed  = zeros(N_S,6)     ;
S_M_bytes  = zeros(N_S,4)     ;
    
S_Wealth_Profile_NB = zeros(M_Simul.p.Max_Age,6,N_S) ; 
S_Wealth_Profile_45 = zeros(M_Simul.p.Max_Age,6,N_S) ; 
S_Wealth_Corr       = zeros(N_S,2)   ; 

% Solve Model
tStart = tic;
M_Simul = Functions_ModelSolution.Aiyagari_Equilibrium(M_Simul);
S_Gamma_timed = toc(tStart);


% Simulate Panel
M_Panel = Functions_MonteCarlo.Simulate_Panel_Timed(M_Simul, M_Panel);


