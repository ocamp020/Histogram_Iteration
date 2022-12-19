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
%    cd()
%    cd './Dropbox/Research/Histogram_Iteration/Julia_Code/Infinite_Horizon_Aiyagari/ '
% Emmanuel's Computer
    % cd()
    % cd 'C:/Users/Emmanuel/Dropbox/RA_Sergio/Histogram_Iteration/Julia_Code/Infinite_Horizon_Aiyagari/ ' % Laptop
    % cd 'D:/Users/Emmanuel/Dropbox/RA_Sergio/Histogram_Iteration/Julia_Code/Infinite_Horizon_Aiyagari/ ' % Desktop
    % cd 'C:/Users/Emmanuel/Dropbox/RA_Sergio/Histogram_Iteration/Julia_Code/Infinite_Horizon_Aiyagari/ '
% Baxter's Computer
    % cd 'D:/Dropbox/Files/Economics-Research/Project-09_SIM/Code/Histogram_Iteration/Julia_Code/Infinite_Horizon_Aiyagari/ '
% Javier's Computer
     %Home
     cd '/Users/cyberdim/Dropbox/WESTERN_ECONOMICS/RA_Baxter/Matlab_Hist_Iter/Infinitely_Lived'
     % Laptop
     %cd '/Volumes/EHDD1/Dropbox/WESTERN_ECONOMICS/RA_Baxter/Matlab_Hist_Iter/Infinitely_Lived'
% Compute Canada Server
%    cd '/scratch/robin370/Histogram_Iteration/Julia_Code/Infinite_Horizon_Aiyagari/ '
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
%rho_zeta: 0.7000
%sigma_zeta: 1.3000
%r: 0.0320
%w: 53.6240
%a_min: 0.0100
%max_iter: 20000
%dist_tol: 5.0000e-06
%dist_tol_Delta: 1.0000e-10
%eta: 0.1000
%Hist_max_iter: 1500
%Hist_tol: 1.0000e-08
%Hist_eta: 0.1000
%c_min: 1.0000e-16

par = Functions_ModelSolution.Parameters(0.94,2.0, ...
            0.963,0.162,0.70,1.30, ...
            0.0320,53.624,0.010,20000, ...
            5e-6,1e-10, ...
            0.10, 1500, 1e-8,0.1, 1e-16);
p=par;
%par = struct('beta',0.94, 'sigma',2.0, 'rho_eps', 0.963, 'sigma_eps',0.162,'rho_zeta',0.70, 'sigma_zeta',1.30,...
%    'r',0.0320, 'w',53.624,...
%    'a_min',0.010,...
%    'max_iter',20000,...
%    'dist_tol',5e-6,'dist_tol_Delta',1e-10,'eta',0.10,...
%    'Hist_max_iter',1500,'Hist_tol',1e-8,'Hist_eta',0.1,...
%    'c_min',1e-16); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model
% p: [1×1 struct]
% a_max: 120000
% theta_a: 4.5000
% theta_a_f: 4.5000
% n_a: 250
% n_a_fine: 500
% n_zeta = 7;
% n_eps = 11;
% method = 1; % 1 for Kronecker and 2 for loops in expectation of PFI 
% read_flag = false; % Boolean for reading results from file. 

Model = Functions_ModelSolution.ModelStructure(par, ...
    120000, 4.5, 4.5, 250, 500, 7,11,1,false);
M=Model;



%Model = struct('p',par, 'a_max', 120000, 'theta_a', 4.5, 'theta_a_f', 4.5, 'n_a',250, 'n_a_fine',500);
%Model.a_grid = VFI_Toolbox.Make_Grid(Model.n_a, Model.theta_a, par.a_min, Model.a_max,"Poly");
%Model.a_grid_fine = VFI_Toolbox.Make_Grid(Model.n_a_fine, Model.theta_a_f, par.a_min, Model.a_max,"Poly");
%Model.n_cut_fine = double(VFI_Toolbox.Grid_Inv(20000,Model.n_a_fine,Model.theta_a_f, par.a_min, Model.a_max,"Poly"));
% Interest rate process
%Model.n_zeta = 7;
%Model.MP_zeta = struct('N',{}, 'grid',{}, 'Pi',{}, 'PDF',{},'CDF',{});
%[Model.MP_zeta(1).N,Model.MP_zeta(1).grid,Model.MP_zeta(1).Pi,Model.MP_zeta(1).PDF,Model.MP_zeta(1).CDF] = VFI_Toolbox.Tauchen86(par.rho_zeta,par.sigma_zeta,Model.n_zeta,1.96);
%Model.zeta_ref = 1/dot(exp(Model.MP_zeta(1).grid),Model.MP_zeta(1).PDF);
%Model.zeta_grid = Model.zeta_ref.*exp(Model.MP_zeta(1).grid);
% Labor productivity process
%Model.n_eps = 11;
%Model.MP_eps = struct('N',{}, 'grid',{}, 'Pi',{}, 'PDF',{},'CDF',{});
%[Model.MP_eps(1).N,Model.MP_eps(1).grid,Model.MP_eps(1).Pi,Model.MP_eps(1).PDF,Model.MP_eps(1).CDF] = VFI_Toolbox.Rouwenhorst95(par.rho_eps,par.sigma_eps,Model.n_eps);
%Model.eps_ref = 1/dot(exp(Model.MP_eps(1).grid),Model.MP_eps(1).PDF);
%Model.eps_grid = Model.eps_ref.*exp(Model.MP_eps(1).grid);
% State matrices
%Model.a_mat = repmat(transpose(Model.a_grid),[1,Model.n_eps,Model.n_zeta]);
%Model.eps_mat = repmat(Model.eps_grid,[Model.n_a,1,Model.n_zeta]);
%Model.zeta_mat = repmat(reshape(Model.zeta_grid,[1,1,Model.n_zeta]),Model.n_a, Model.n_eps,1);
%Model.a_mat_fine =  repmat(transpose(Model.a_grid_fine),[1,Model.n_eps,Model.n_zeta]);
%Model.eps_mat_fine = repmat(Model.eps_grid,[Model.n_a_fine,1,Model.n_zeta]);
%Model.zeta_mat_fine = repmat(reshape(Model.zeta_grid,[1,1,Model.n_zeta]),Model.n_a_fine, Model.n_eps,1);
% Value and policy functions 
%Model.V = zeros([Model.n_a, Model.n_eps, Model.n_zeta]); % Value function
%Model.G_ap = zeros([Model.n_a, Model.n_eps, Model.n_zeta]); % Policy Function for Capital/Assets
%Model.G_c = zeros([Model.n_a, Model.n_eps, Model.n_zeta]); % Policy Function for Consumption
%Model.V_fine = zeros([Model.n_a_fine, Model.n_eps, Model.n_zeta]); % Value function
%Model.G_ap_fine = zeros([Model.n_a_fine, Model.n_eps, Model.n_zeta]); % Policy Function for Capital/Assets
%Model.G_c_fine = zeros([Model.n_a_fine, Model.n_eps, Model.n_zeta]); % Policy Function for Consumption
% Distribution
%Model.Gamma = (1/(Model.n_cut_fine*Model.n_eps*Model.n_zeta)).*cat(1,ones(Model.n_cut_fine,Model.n_eps,Model.n_zeta), zeros(Model.n_a_fine-Model.n_cut_fine, Model.n_eps, Model.n_zeta));
%Model.H_ind = rand(Model.n_a_fine, Model.n_eps, Model.n_zeta); % Index for discretization of savings choice  
%Model.H_omega_lo = rand(Model.n_a_fine, Model.n_eps, Model.n_zeta, Model.n_eps, Model.n_zeta); % Index for discretization of savings choice  
%Model.H_omega_hi = rand(Model.n_a_fine, Model.n_eps, Model.n_zeta, Model.n_eps, Model.n_zeta);
% Misc 
%Model.method = 1; % 1 for Kronecker and 2 for loops in expectation of PFI 
%Model.read_flag = false; % Boolean for reading results from file. 
%M = Model;
% End of model


% Load functions in Functions_ModelSolution (solve the model and find stationary distribution)
% Load functions in Functions_Montecarlo (simulate panels of individual agents)

% Execute model solution
disp('===============================================')
disp('Solving Aiyagari with EGM-Histogram(loop)')
%format longG
%
%
% 
tic
M_Aiyagari = Functions_ModelSolution.Aiyagari_Equilibrium(Model);
toc

% # Get stats and graphs for the solution of the model 
%run('CalculateMoments_Solution.m')


disp('===============================================')



% Run Draft Moments for Graphs and Tables 
run("Draft_Results.m")



disp('=============================================== End of Script ===============================================')








