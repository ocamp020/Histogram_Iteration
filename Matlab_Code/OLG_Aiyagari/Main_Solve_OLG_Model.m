% % Computing Longitudinal Moments for Heterogeneous Agent Models
% Sergio Ocampo, Baxter Robinson, Emmanuel Murray Leclair and Javier Martinez
% April 2022
% Aiyagari economy: 
%       1. Infinitely lived agents
%       2. Inelastic labor supply
%       3. Stochastic rate of returns
% This scripts computes longitudinal moments for the model
% 
% Solve:   V(ζ,ϵ,a) = max{ ((1+r(ζ))a+wϵ̄ϵ-a')^(1-σ)/(1-σ) +beta*E[V(ζ',ϵ',a')|ζ,ϵ] }
%           log(ϵ') = ρ_ϵ*log(ϵ) + η_ϵ; η_ϵ~N(0,σ_ϵ);
%           r(ζ)    = exp(ζ)r⋆    
%           log(ζ') = ρ_ζ*log(ζ) + η_ζ; η_ζ~N(0,σ_ζ); 
% The constant ϵ̄ guarantees that E[ϵ]=1 and so aggregate labor L=E[ϵ]=1

%% Change to your home directory 
% Sergio's Computer 
  % cd()
   %cd("./Dropbox/Research/Histogram_Iteration/Julia_Code/OLG_Aiyagari/")
% Emmanuel's Computer
    % cd()
    % cd("C:/Users/Emmanuel/Dropbox/RA_Sergio/Histogram_Iteration/Julia_Code/OLG_Aiyagari/") # Laptop
    % cd("D:/Users/Emmanuel/Dropbox/RA_Sergio/Histogram_Iteration/Julia_Code/OLG_Aiyagari/") # Desktop
    % cd("C:/Users/Emmanuel/Dropbox/RA_Sergio/Histogram_Iteration/Julia_Code/OLG_Aiyagari/")
% Baxter's Computer
    %  cd("D:/Dropbox/Files/Economics-Research/Project-09_SIM/Code/Histogram_Iteration/Julia_Code/OLG_Aiyagari/")
% Javier's Computer
     %Home
     cd '/Users/cyberdim/Dropbox/WESTERN_ECONOMICS/RA_Baxter/Matlab_Hist_Iter/OLG'
     % Laptop
     %cd '/Volumes/EHDD1/Dropbox/WESTERN_ECONOMICS/RA_Baxter/Matlab_Hist_Iter/OLG'
% Compute Canada Server
%    cd '/scratch/robin370/Histogram_Iteration/Julia_Code/OLG_Aiyagari/ '
clc
clear


% Make auxiliary directories
mkdir Fig_Folder
mkdir File_Folder
mkdir Hist_Folder
mkdir MC_Folder


disp(" ")
disp("------------------------")
disp("Aiyagari in Matlab")
disp("PWD: ")
disp(pwd)
disp("Solve Aiyagari model with EGM and the histogram method")
disp("------------------------")
disp(" ")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters

%beta: 0.9400
%sigma: 2
%rho_eps: 0.9630
%sigma_eps: 0.1620
%r: 0.0379
%w: 53.6240
%a_min: 1.0000e-03
%Hist_max_iter: 1000
%Hist_tol: 1.0000e-06
%Hist_eta: 0
%c_min: 1.0000e-16
%Max_Age: 81
%Surv_Pr: [81×1 double]
%Age_Pi: [81×81 double]
%Age_PDF: [81×1 double]

par = Functions_ModelSolution.Parameters(0.94, 2.0, 0.963,0.162,...
    0.0379, 53.624,1e-3,1000,1e-6,0.00, 1e-16,81);
p=par;

%par = struct('beta',0.94, 'sigma',2.0, 'rho_eps', 0.963, 'sigma_eps',0.162,...
%    'r',0.0379, 'w',53.624,...
%    'a_min',1e-3,...
%'Hist_max_iter',1000,'Hist_tol',1e-6,'Hist_eta',0.00,...
%    'c_min',1e-16, ...
%    'Max_Age',81);
%par.Surv_Pr = Setup_Demographics.Survival_Probabilities_Bell_Miller(par.Max_Age);
%par.Age_Pi = Setup_Demographics.Age_Transition(par.Max_Age, par.Surv_Pr);
%par.Age_PDF = Setup_Demographics.Age_Distribution(par.Age_Pi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model

%par: [1×1 struct]
%a_max: 10000
%theta_a: 3.5000
%theta_a_f: 3.5000
%n_a: 250
%n_a_fine: 500
%n_eps: 11
%read_flag: 0

Model = Functions_ModelSolution.ModelStructure(par,10000,3.5,3.5,250, 500,11,0 );
M= Model;
%Model = struct('p',par, ...
%    'a_max', 10000, 'theta_a', 3.5, 'theta_a_f', 3.5, 'n_a',250, 'n_a_fine',500);
%Model.a_grid = VFI_Toolbox.Make_Grid(Model.n_a, Model.theta_a, par.a_min, Model.a_max,"Poly");
%Model.a_grid_fine = VFI_Toolbox.Make_Grid(Model.n_a_fine, Model.theta_a_f, par.a_min, Model.a_max,"Poly");

% Labor productivity process
%Model.n_eps = 11; % Size of ϵ_grid
%Model.MP_eps = struct('N',{}, 'grid',{}, 'Pi',{}, 'PDF',{},'CDF',{});
%[Model.MP_eps(1).N,Model.MP_eps(1).grid,Model.MP_eps(1).Pi,Model.MP_eps(1).PDF,Model.MP_eps(1).CDF] = VFI_Toolbox.Rouwenhorst95(par.rho_eps,par.sigma_eps,Model.n_eps); % Markov Process for epsilon
%Model.eps_ref = 1.038479216975849/dot(exp(Model.MP_eps(1).grid),Model.MP_eps(1).PDF); % Reference level for labor efficiency 
%Model.eps_grid = Model.eps_ref.*exp(Model.MP_eps(1).grid); % Grid in levels

% Labor productivity process - Life Cycle 
%Model.age_vec = transpose(1:1:par.Max_Age) ;
%Model.log_xi_grid = (60*(Model.age_vec-1) - ((Model.age_vec-1).^2))/1800 ; % Process peaks at age 50, and by age 80 gives the same income as when newborn
%Model.xi_ref = 1/sum(exp(Model.log_xi_grid).*par.Age_PDF) ; % Reference level for labor efficiency 
%Model.xi_grid = Model.xi_ref*exp(Model.log_xi_grid) ; %Grid in levels

% State matrices
%Model.a_mat = repmat(transpose(Model.a_grid),[1,Model.n_eps,par.Max_Age]);
%Model.eps_mat = repmat(Model.eps_grid,[Model.n_a,1,par.Max_Age]);
%Model.xi_mat = repmat(reshape(Model.xi_grid,[1,1,par.Max_Age]),Model.n_a, Model.n_eps,1);
%Model.a_mat_fine =  repmat(transpose(Model.a_grid_fine),[1,Model.n_eps,par.Max_Age]);
%Model.eps_mat_fine = repmat(Model.eps_grid,[Model.n_a_fine,1,par.Max_Age]);
%Model.xi_mat_fine = repmat(reshape(Model.xi_grid,[1,1,par.Max_Age]),Model.n_a_fine, Model.n_eps,1);
%Model.a_mat_aeps = repmat(transpose(Model.a_grid),1,Model.n_eps);

% Labor income matrices
%Model.y_mat = par.w*Model.eps_mat.*Model.xi_mat ; 
%Model.y_mat_fine = par.w*Model.eps_mat_fine.*Model.xi_mat_fine ; 

% Value and Policy Functions
%Model.V = zeros([Model.n_a, Model.n_eps, par.Max_Age]); % Value function
%Model.G_ap = zeros([Model.n_a, Model.n_eps, par.Max_Age]); % Policy Function for Capital/Assets
%Model.G_c = zeros([Model.n_a, Model.n_eps, par.Max_Age]); % Policy Function for Consumption
%Model.V_fine = zeros([Model.n_a_fine, Model.n_eps, par.Max_Age]); % Value function
%Model.G_ap_fine = zeros([Model.n_a_fine, Model.n_eps, par.Max_Age]); % Policy Function for Capital/Assets
%Model.G_c_fine = zeros([Model.n_a_fine, Model.n_eps, par.Max_Age]); % Policy Function for Consumption

% Distribution Guess
%Model.n_cut_fine = double(VFI_Toolbox.Grid_Inv(1000,Model.n_a_fine,Model.theta_a_f, par.a_min, Model.a_max,"Poly")); %  Index just below 1000
%Model.Gamma_a_guess = (1/((Model.n_cut_fine))).*([ones([Model.n_cut_fine, Model.n_eps, par.Max_Age]) ; zeros([Model.n_a_fine - Model.n_cut_fine, Model.n_eps, par.Max_Age])]);
%Model.Gamma = Model.Gamma_a_guess.*repmat(Model.MP_eps(1).PDF, Model.n_a_fine, 1, par.Max_Age).*repmat(reshape(par.Age_PDF,[1,1,par.Max_Age]),Model.n_a_fine, Model.n_eps,1);

% Matrices for discretization of policy functions
%Model.H_ind = zeros([Model.n_a_fine, Model.n_eps, par.Max_Age]); % Index for discretization of savings choice  
%Model.H_omega_lo_s = zeros(Model.n_a_fine, Model.n_eps, par.Max_Age, Model.n_eps); 
%Model.H_omega_hi_s = zeros(Model.n_a_fine, Model.n_eps, par.Max_Age, Model.n_eps);
%Model.H_omega_lo_d = zeros(Model.n_a_fine, Model.n_eps, par.Max_Age, Model.n_eps);
%Model.H_omega_hi_d = zeros(Model.n_a_fine, Model.n_eps, par.Max_Age, Model.n_eps);
% Misc 
%Model.read_flag = false; % Boolean for reading results from file. 

% End of model


% Execute model solution
disp('===============================================')
disp('Solving Aiyagari with EGM-Histogram(loop)')
%format longG
%
%
% 
%

tStart = tic;
M_Aiyagari = Functions_ModelSolution.Aiyagari_Equilibrium(Model);
tEnd = toc(tStart);

% Get stats and graphs for the solution of the model 
run('CalculateMoments_Solution.m')

disp('===============================================')




% Run Draft Moments for Graphs and Tables 
% include("Draft_Results.jl")



disp('=============================================== End of Script ===============================================')








