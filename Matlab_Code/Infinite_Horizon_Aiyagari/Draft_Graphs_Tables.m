%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Graphs and Tables for Draft
%% Loads results generated by Draft_Results



% Set up model structures
Model1 = Functions_ModelSolution.ModelStructure(par, ...
    120000, 4.5, 4.5, 250, 500, 7,11,1,false);
M1=Model1;
clear Model1
M_Simul = Functions_ModelSolution.model2modify(M1, 'read_flag',1 );
M_Panel = Functions_MonteCarlo.PanelStructure(500000, 10, 1000, 1000, 1e-2, 3489398);
M_Panel = Functions_MonteCarlo.panelmodel2modify(M_Panel,'N_Panel',1000000);


% Load results from csv files
H_grid_size = [250 250 500 1000 5000];
n_H = length(H_grid_size);
pct_list = [90;95;99;99.9;99.99];
S_sample = [50000 250000 500000 1000000 10000000];
N_S = length(S_sample);

H_Gamma_timed = csvread(strcat(Hist_Folder, '/H_G_timed.csv'));
H_Gamma_bytes = csvread(strcat(Hist_Folder, '/H_G_bytes.csv'));
H_M_timed = reshape(csvread(strcat(Hist_Folder, '/H_M_timed.csv')), n_H, 3);
H_M_bytes = reshape(csvread(strcat(Hist_Folder, '/H_M_bytes.csv')), n_H, 3);

H_a_grid = reshape(csvread(strcat(Hist_Folder, '/H_a_grid.csv')), H_grid_size(end), n_H);
H_Gamma_a = reshape(csvread(strcat(Hist_Folder, '/H_G_a.csv')), H_grid_size(end), n_H);

H_Wealth_Stats = reshape(csvread(strcat(Hist_Folder, '/H_Wealth_Stats.csv')), n_H, 6);
H_Wealth_Share = reshape(csvread(strcat(Hist_Folder, '/H_Wealth_Share.csv')), n_H, 5);
H_Pareto_Coeff = csvread(strcat(Hist_Folder, '/H_Pareto_Coeff.csv'));
H_Decile = reshape(csvread(strcat(Hist_Folder, '/H_Decile.csv')), 11, 2, n_H);
H_Decile_Tr = reshape(csvread(strcat(Hist_Folder, '/H_Decile_Tr.csv')), 10, 10, n_H);
H_Cons_Corr = csvread(strcat(Hist_Folder, '/H_Cons_Corr.csv'));
H_A_Corr = csvread(strcat(Hist_Folder, '/H_A_Corr.csv'));
H_eps_Corr = csvread(strcat(Hist_Folder, '/H_eps_Corr.csv'));
H_zeta_Corr = csvread(strcat(Hist_Folder, '/H_z_Corr.csv'));

S_M_timed = reshape(csvread(strcat(MC_Folder, '/S_M_timed.csv')), N_S, 4);
S_M_bytes = reshape(csvread(strcat(MC_Folder, '/S_M_bytes.csv')), N_S, 4);
S_Gamma_timed = csvread(strcat(MC_Folder, '/S_G_timed.csv'));


S_Wealth_Sample_0 = reshape(csvread(strcat(MC_Folder, '/S_Wealth_Sample.csv')), N_S-1, S_sample(N_S-1));
S_Wealth_Sample_Large = reshape(csvread(strcat(MC_Folder, '/S_Wealth_Sample_Large.csv')), 1, S_sample(N_S));
S_Wealth_Sample = zeros(N_S, S_sample(N_S));
S_Wealth_Sample(1:N_S-1, 1:S_sample(N_S-1)) = S_Wealth_Sample_0;
S_Wealth_Sample(N_S, :) = S_Wealth_Sample_Large;


S_Wealth_Stats = reshape(csvread(strcat(MC_Folder, '/S_Wealth_Stats.csv')), N_S, 6);
S_Wealth_Share = reshape(csvread(strcat(MC_Folder, '/S_Wealth_Share.csv')), N_S, 5);
S_Pareto_Coeff = csvread(strcat(MC_Folder, '/S_Pareto_Coeff.csv'));
S_Decile = reshape(csvread(strcat(MC_Folder, '/S_Decile.csv')), 11, N_S);
S_Decile_Tr = reshape(csvread(strcat(MC_Folder, '/S_Decile_Tr.csv')), 10, 10, N_S);
S_Cons_Corr = csvread(strcat(MC_Folder, '/S_Cons_Corr.csv'));
S_A_Corr = csvread(strcat(MC_Folder, '/S_A_Corr.csv'));
S_eps_Corr = csvread(strcat(MC_Folder, '/S_eps_Corr.csv'));
S_zeta_Corr = csvread(strcat(MC_Folder, '/S_z_Corr.csv'));


Mat_Top_Stats = {'Top Wealth Shares + Pareto Tail', '', '', '', '', ''; ...
                'Histogram', '', '', '', '', ''; ...
                'Grid Size', H_grid_size'; ...
                'Top 0.1%', H_Wealth_Share(:,4)'; ...
                'Top 1%', H_Wealth_Share(:,3)'; ...
                'Top 10%', H_Wealth_Share(:,1)'; ...
                'Pareto', H_Pareto_Coeff'; ...
                'Av. Assets', H_Wealth_Stats(:,end)'; ...
                '-', '-', '-', '-', '-', '-'; ...
                'Total Time', H_Gamma_timed' + H_M_timed(:,1)'; ...
                'Model Time', H_Gamma_timed'; ...
                'Moment Time', H_M_timed(:,1)'; ...
                '-', '-', '-', '-', '-', '-'; ...
                '-', '-', '-', '-', '-', '-'; ...
                '-', '-', '-', '-', '-', '-'; ...
                'Simulation', '', '', '', '', ''; ...
                'Grid Size', '50k', '250k', '500k', '1M', '10M'; ...
                'Top 0.1%', S_Wealth_Share(:,4)'; ...
                'Top 1%', S_Wealth_Share(:,3)'; ...
                'Top 10%', S_Wealth_Share(:,1)'; ...
                'Pareto', S_Pareto_Coeff'; ...
                'Av. Assets', S_Wealth_Stats(:,end)'; ...
                '-', '-', '-', '-', '-', '-'; ...
                'Total Time',   S_Gamma_timed + S_M_timed(:,1)' + S_M_timed(:,2)'  ; ...
                'Model Time',    repeat((S_Gamma_timed(1)),1,N_S)   ; ...
                'Simul Time',    S_M_timed(:,1)'             ; ...
                'Moment Time',   S_M_timed(:,2)'             ; ...
                '-', '-', '-', '-', '-', '-';};

writematrix(Mat_Top_Stats,'./Fig_Folder/Table_Top_Stats.csv')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up colors
%(LINE 132 IN JULIA)