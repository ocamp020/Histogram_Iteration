
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions for Monte-Carlo simulation 


    classdef Functions_MonteCarlo
        methods(Static)

            % Adapted from Inverse transform sampling for Panel Simulation
            % pdfrnd(x, px, sampleSize): return a random sample of size sampleSize from 
            % the pdf px defined on the domain x. Taken from 
            % Joshua Stough (2022). 
            % Random Sample from Discrete PDF 
            % https://www.mathworks.com/matlabcentral/fileexchange/37698-random-sample-from-discrete-pdf
            % MATLAB Central File Exchange. Retrieved August 16, 2022. 

            function [X] = pdfrnd(px, sampleSize)
                rnd = rand(sampleSize, 1); % Draw a random number from Uniform distribution [0,1]
                % Find the cdf given our pdf (bc we are looking for the
                % cdf^{-1}
                if isrow(px)
                    px = px';
                else
                    px = px;
                end
                px = px/sum(px);
                cdf = cumsum(px);
                % Inversion 
                X = zeros(sampleSize,1);
                for k=1:sampleSize
                    I = find(cdf>=rnd(k),1,'first');
                    X(k)=I;
                end
            end % End of pdfrnd



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            function  Model_Panel = Model_Panel(N_Panel, T_Panel, T_Simul, N_Min, Simul_tol, rng_seed )
                
                Model_Panel = struct('N_Panel', N_Panel , 'T_Panel', T_Panel, 'T_Simul', T_Simul, 'N_Min', N_Min, 'Simul_tol', Simul_tol, 'rng_seed', rng_seed); 
                %Model_Panel = struct('N_Panel', 500000 , 'T_Panel', 10, 'T_Simul', 1500, 'N_Min', 1000, 'Simul_tol', 1e-2, 'rng_seed', 297835398); 
                % Panel Output
                Model_Panel.a_mat    = zeros([Model_Panel.N_Panel, Model_Panel.T_Panel]);
                Model_Panel.c_mat    =  zeros([Model_Panel.N_Panel, Model_Panel.T_Panel]);
                Model_Panel.eps_mat  = zeros([Model_Panel.N_Panel, Model_Panel.T_Panel]);
                Model_Panel.h_mat    =  zeros([Model_Panel.N_Panel, Model_Panel.T_Panel]);
                Model_Panel.t_mat    = zeros([Model_Panel.N_Panel,1]);
            
            end % End of function Model_Panel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            function updated_model_structure = model_panel2modify(model_structure, parameter_requested, parameter_new_value)
                % This function checks if the requested parameter is in the
                % model structure. If so, requests the function 'ModelStructure' for 
                % an updated model structure. If the requested parameter is not
                % there, displays an error sign.
        
                % This is what I though to use to overcome the @kw macro/structure
                % in Julia
        
                % The requested parameter should be written as: 'parameter_requested' .
                
                old = model_structure;
                if any(strcmp(fieldnames(model_structure),parameter_requested))
                    old.(sprintf(parameter_requested)) = parameter_new_value;
                    values = struct2cell(old);
                    %fields_inner
                    updated_model_structure = Functions_ModelSolution.ModelStructure(values{1:6});
                else
                    error('Parameter not in the structure')
                end
                clear old
                return
                
            end % End of function model_panel2modify


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            function Model_Cohort = Model_Cohort(N_Panel, T_Panel, rng_seed )
            
                Model_Cohort = struct('N_Panel', N_Panel, 'T_Panel', T_Panel, 'rng_seed', rng_seed);
                %Model_Cohort = struct('N_Panel', 500000, 'T_Panel', 100, 'rng_seed', 297835398);
                Model_Cohort.a_mat    = zeros([Model_Cohort.N_Panel, Model_Cohort.T_Panel]);
                Model_Cohort.c_mat    =  zeros([Model_Cohort.N_Panel, Model_Cohort.T_Panel]);
                Model_Cohort.eps_mat  = zeros([Model_Cohort.N_Panel, Model_Cohort.T_Panel]);
                Model_Cohort.h_mat    =  zeros([Model_Cohort.N_Panel, Model_Cohort.T_Panel]);
                Model_Cohort.t_mat    = zeros([Model_Cohort.N_Panel,1]);
            
            
            end % End of function Model_Cohort
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            function updated_model_structure = model_cohort2modify(model_structure, parameter_requested, parameter_new_value)
                            % This function checks if the requested parameter is in the
                            % model structure. If so, requests the function 'ModelStructure' for 
                            % an updated model structure. If the requested parameter is not
                            % there, displays an error sign.
                    
                            % This is what I though to use to overcome the @kw macro/structure
                            % in Julia
                    
                            % The requested parameter should be written as: 'parameter_requested' .
                            
                            old = model_structure;
                            if any(strcmp(fieldnames(model_structure),parameter_requested))
                                old.(sprintf(parameter_requested)) = parameter_new_value;
                                values = struct2cell(old);
                                %fields_inner
                                updated_model_structure = Functions_ModelSolution.ModelStructure(values{1:3});
                            else
                                error('Parameter not in the structure')
                            end
                            clear old
                            return
                
            end % End of function model_cohort2modify


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Function Simulate_Panel
            function ModelStructure_Output = Simulate_Panel(ModelStructure,PanelModel_Structure) % 
                
                % ModelStructure is the input structure (M in Julia code)
                % PanelModel_Structure is the structure of the panel parameters (M_P in Julia code)
                
                % ----------- Output variables ----------- 
                % Output1 = PM (updated version) 
                % Output1 is a Panel Model Structure 
                % ----------------------------------------  


                % Local Variables
                M = ModelStructure;
                PM = PanelModel_Structure;

                p = M.p;
                MP_eps = M.MP_eps;
                a_grid = M.a_grid;
                a_grid_fine = M.a_grid_fine;
                a_max = M.a_max;
                n_eps = M.n_eps;
                Gamma = M.Gamma;
                G_ap_fine = M.G_ap_fine;
                G_c_fine = M.G_c_fine;
                n_a_fine = M.n_a_fine;
                theta_a_f = M.theta_a_f;
                eps_grid = M.eps_grid;

                a_min = p.a_min;
                Surv_Pr = p.Surv_Pr;
                Max_Age = p.Max_Age;
                

                N_Panel = PM.N_Panel;
                T_Panel = PM.T_Panel;
                T_Simul = PM.T_Simul;
                

                %Seed_Flag = 1;
                
                % Initialize Seed
                rng(PM.rng_seed);

                % Initialize cross-section vectors
                a_vec = zeros(N_Panel,1);
                c_vec = zeros(N_Panel,1);
                eps_vec = zeros(N_Panel,1);
                h_vec = zeros(N_Panel,1);

                % PDFs for epsilon and zeta
                Gamma_eps  = MP_eps.Pi;
                
                % Set median epsilon for newborns
                med_eps = int64(round(n_eps/2));

                % Set initial wealth to (close to) $1k
                b_ind = transpose(find(a_grid_fine>=1));
                b_ind = b_ind(1);
                b = a_grid_fine(b_ind);

                % Censor savings
                G_ap_fine = max(min(G_ap_fine,a_max), a_min);

                % Draw initial conditions from stationary distribution (cross-section)
                fprintf(' Initializing simulation \n')

                % Draw eps(i) and zeta(i) from (stationary) marginal
                % distribution
                eps_pdf = transpose(squeeze(sum(Gamma, [1,3])));
                eps_vec = Functions_MonteCarlo.pdfrnd(eps_pdf, N_Panel);
                h_pdf = squeeze(sum(M.Gamma, [1,2]));
                h_vec = Functions_MonteCarlo.pdfrnd(h_pdf, N_Panel);
                clear eps_pdf h_pdf

                % Draw a(i) from conditional distribution 

               for i=1:N_Panel
                   if h_vec(i)==1
                       a_vec(i)=b;
                   else
                       Gamma_aux = Gamma(:,eps_vec(i), h_vec(i))/sum(Gamma(:,eps_vec(i), h_vec(i)));
                       a_vec(i) = a_grid_fine(Functions_MonteCarlo.pdfrnd(Gamma_aux, 1));
                   end
                end
                
                
                %*************************************************************
                %*************************************************************
                % Importing julia's a_vec, eps_vec, zeta_vec 
                % bc they are random vectors so I just want to check the
                % code has the same output given that it 
                % works with the same input numbers
                %eps_vec = readmatrix('/Users/cyberdim/Dropbox/WESTERN_ECONOMICS/RA_Baxter/Matlab_Hist_Iter/OLG/File_Folder/eps_vec_Julia.csv');
                %h_vec = readmatrix('/Users/cyberdim/Dropbox/WESTERN_ECONOMICS/RA_Baxter/Matlab_Hist_Iter/OLG/File_Folder/h_vec_Julia.csv');
                %a_vec = readmatrix('/Users/cyberdim/Dropbox/WESTERN_ECONOMICS/RA_Baxter/Matlab_Hist_Iter/OLG/File_Folder/a_vec_Julia.csv');
                %*************************************************************
                %*************************************************************

                fprintf('E[a_0] = $%8.2f k , E[eps_0]= %8.3f , E[h_0] = %8.3f, max(a_0) = %8.2f \n', ...
                    round(mean(a_vec),2), round(mean(eps_grid(eps_vec)),3), ...
                    round(mean(h_vec),1), round(max(a_vec),2)  )

                

                % Iterate Forward
                disp('Iterating Panel')
                for t=1:T_Simul
                    if mod(t,50)==0
                        fprintf(' Simulation Period %i \n', t)
                        fprintf(' E[a_t] = %4.2fk E[eps_t]=%3.3f E[h_t]=%4.1f, max(a_t) = %8.2f  \n', round(mean(a_vec),2), round(mean(eps_grid(eps_vec)),3), ...
                    round(mean(h_vec),1), round(max(a_vec),2) )
                    end
                    % Simulate each dinasty
                    for i=1:N_Panel
                        Death_Draw = Functions_MonteCarlo.pdfrnd([Surv_Pr(h_vec(i)), 1-Surv_Pr(h_vec(i)) ],1);

                        if Death_Draw ==2 % Agent dies
                            a_vec(i)=b;
                            eps_vec(i)=med_eps;
                            h_vec(i)=1;
                        else % Agent lives
                            % Compute future assets: Linear interpolation (manually)
                            if (a_min<a_vec(i)) && (a_vec(i)<a_max)
                                i_lo = min(n_a_fine-1, VFI_Toolbox.Grid_Inv(a_vec(i), n_a_fine, theta_a_f, a_min, a_max, "Poly" ));
                                omega = min(1, max(0,(a_vec(i) - a_grid_fine(i_lo))/(a_grid_fine(i_lo+1) - a_grid_fine(i_lo))))  ;
                                a_vec(i) = omega*G_ap_fine(i_lo, eps_vec(i), h_vec(i)) + (1-omega)*G_ap_fine(i_lo+1, eps_vec(i), h_vec(i))  ;
                            elseif a_vec(i)==a_max
                                 a_vec(i) = G_ap_fine(end, eps_vec(i), h_vec(i)) ;
                            elseif a_vec(i)==a_min
                                a_vec(i) = G_ap_fine(1, eps_vec(i), h_vec(i)) ;
                            else
                            error('Error in simulation: assets not working')
                            end
                            a_vec(i) = min(a_max, max(a_min,a_vec(i)));
                            % Compute future eps(i)
                            eps_pdf = transpose(Gamma_eps(eps_vec(i),:));
                            eps_vec(i) = Functions_MonteCarlo.pdfrnd(eps_pdf,1);
                            clear eps_pdf
                            % Update age
                            h_vec(i) = h_vec(i)+1;
                        end
                     %Compute consumption only for relevant periods 
                        if t>=T_Simul - (T_Panel - 1 )
                            if (a_min < a_vec(i)) && (a_vec(i) < a_max)
                                i_lo = min(n_a_fine-1, VFI_Toolbox.Grid_Inv(a_vec(i), n_a_fine, theta_a_f, a_min, a_max, "Poly" ));
                                omega = min(1, max(0,(a_vec(i) - a_grid_fine(i_lo))/(a_grid_fine(i_lo+1) - a_grid_fine(i_lo)   )))  ;
                                c_vec(i) = omega*G_c_fine(i_lo, eps_vec(i), h_vec(i)) + (1-omega)*G_c_fine(i_lo+1, eps_vec(i), h_vec(i) )  ;
                            elseif a_min==a_vec(i)
                                c_vec(i) = G_c_fine(1,eps_vec(i), h_vec(i));
                            else
                                c_vec(i) = G_c_fine(end,eps_vec(i), h_vec(i));
                            end
                        end
                    end
                    
                  % Save results in panel
                        if t>= T_Simul - (T_Panel - 1)
                            PM.a_mat(:, t-(T_Simul-T_Panel)) = a_vec ; % Assets
                            PM.eps_mat(:, t-(T_Simul-T_Panel)) = eps_vec ; % Assets
                            PM.h_mat(:, t-(T_Simul-T_Panel)) = h_vec ; % Assets
                            PM.c_mat(:, t-(T_Simul-T_Panel)) = c_vec ; % Assets
                        end
                end
                

                ModelStructure_Output = PM;
                return
               
            end % End of Function Simulate_Panel





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Simulate Cohort

            function Updated_PanelModel_Structure = Simulate_Cohort(ModelStructure, PanelModel_Structure_Cohort, age_0) 
                
                
                % c_flag is a Boolean variable

                % ----------- Output contains 4 variables ----------- 
                % Output1 = a_mat
                % Output2 = c_mat
                % Output3 = eps_mat
                % Output4 = h_mat
                % ----------------------------------------  
                
                % Local Variables
                M = ModelStructure;
                M_P_cohort = PanelModel_Structure_Cohort;
                 p = M.p;

                MP_eps = M.MP_eps;
                a_grid = M.a_grid;
                eps_grid = M.eps_grid;
                a_grid_fine = M.a_grid_fine;
                a_max = M.a_max;
                n_eps = M.n_eps;
                Gamma = M.Gamma;
                G_ap_fine = M.G_ap_fine;
                G_c_fine = M.G_c_fine;
                n_a_fine = M.n_a_fine;
                theta_a_f = M.theta_a_f;
                
                a_min = p.a_min;
                Surv_Pr = p.Surv_Pr;
                Max_Age = p.Max_Age;


                N_Panel = M_P_cohort.N_Panel;

                % Max Cohort Age
                Max_C_Age = Max_Age - (age_0-1);

                % Initialize seed
                rng(M_P_cohort.rng_seed);

                % Initialize cross-section vector
                
                a_mat    = zeros([N_Panel, Max_C_Age]);
                c_mat    =  zeros([N_Panel, Max_C_Age]);
                eps_mat  = zeros([N_Panel, Max_C_Age]);
                h_mat    = zeros([N_Panel,Max_C_Age]);

                % PDF for epsilon
                Gamma_eps = MP_eps.Pi;

                % Set median epsilon for newborns
                med_eps = int64(round(n_eps/2));

                % Set initial wealth to (close to) $1k
                b_ind = transpose(find(a_grid_fine>=1));
                b_ind = b_ind(1);
                b = a_grid_fine(b_ind);

                % Censor savings
                G_ap_fine = max(min(G_ap_fine, a_max), a_min);
                
                % Draw initial conditions from stationary distribution (cross-section)
                disp("Initializing Simulation")

                % Draw epsilon(i) from (stationary) marignal distribution
                pdf_epsmat = transpose(sum(Gamma(:,:,age_0)/sum(sum(sum(Gamma(:,:,age_0)))),[1,3]));
                eps_mat(:,1) = Functions_MonteCarlo.pdfrnd(pdf_epsmat,N_Panel);
                h_mat(:,1)=1;
                % Draw a(i) from conditional distribution
                if age_0==1
                    a_mat(:,1) = b;
                    c_mat(:,1) = G_c_fine(b_ind,med_eps,1);
                else
                    for i=1:N_Panel
                        Gamma_aux = Gamma(:,eps_mat(i,1),age_0)/sum(Gamma(:,eps_mat(i,1), age_0));
                        a_ind = Functions_MonteCarlo.pdfrnd(Gamma_aux,1);
                        a_mat(i,1) = a_grid_fine(a_ind);
                        c_mat(i,1) = G_c_fine(a_ind, eps_mat(i,1), age_0);
                    end
                end
                v1=round(mean(a_mat(:,1)),2);
                v2=round(mean(eps_grid(eps_mat(:,1))),3);
                v3=round(mean(h_mat(:,1)) ,3);
                v4=round(max(a_mat(:,1)),2);

                fprintf('E[a_0]=%8.2fk  E[eps_0]=%8.3f   E[h_0]=%4.1f  max[a_0]=%8.2f \n', ...
                    v1,v2,v3,v4 )
                clear v1 v2 v3 v4

                % Iterate forward
                disp('   Iterating panel')
                for t=2:Max_C_Age
                    if mod(t,5)==0
                        ind = h_mat(:,t-1)>0; % Select only alive agents
                        fprintf('Simulation Period %i \n',t)
                        fprintf(' E[a_t] = %5.2fk E[eps_t]=%3.3f E[h_t]=%4.1f, max(a_t) = %8.2f  \n', round(mean(a_mat(ind,t-1)),2), round(mean(eps_grid(eps_mat(ind,t-1))),3), ...
                    round(mean(h_mat(:,t-1)),1), round(max(a_mat(ind,t-1)),2))
                    end
                        % Set current age of agents in cohort
                        age = t+(age_0-1);
                        % Simulate each dinasty
                    for i=1:N_Panel
                        % Check if alive // if not just don't update and leave entries as zeros
                        if h_mat(i,t-1)==1
                            Death_Draw = Functions_MonteCarlo.pdfrnd([Surv_Pr(age), 1-Surv_Pr(age) ],1);
                            if Death_Draw ==2 % Agent dies, no need to update 
                                h_mat(i,t)=0;
                            else % Agent lives
                            % Compute future assets: Linear interpolation (manually)
                                if (a_min<a_mat(i,t-1)) && (a_mat(i,t-1)<a_max)
                                    i_lo = min(n_a_fine-1, VFI_Toolbox.Grid_Inv(a_mat(i,t-1), n_a_fine, theta_a_f, a_min, a_max, "Poly" ));
                                    omega = min(1, max(0,(a_mat(i,t-1) - a_grid_fine(i_lo))/(a_grid_fine(i_lo+1) - a_grid_fine(i_lo)   )))  ;
                                    a_mat(i,t) = omega*G_ap_fine(i_lo, eps_mat(i,t-1), age-1) + (1-omega)*G_ap_fine(i_lo+1, eps_mat(i,t-1), age-1 )  ;
                                elseif a_mat(i,t-1)==a_max
                                     a_mat(i,t) = G_ap_fine(end, eps_mat(i,t-1), age-1) ;
                                elseif a_mat(i,t-1)==a_min
                                    a_mat(i,t) = G_ap_fine(1, eps_mat(i,t-1),age-1) ;
                                else
                                error('Error in simulation: assets not working')
                                end
                            a_mat(i,t) = min(a_max, max(a_min,a_mat(i,t)));

                            % Compute future eps(i)
                            eps_mat(i,t) = Functions_MonteCarlo.pdfrnd(transpose(Gamma_eps(eps_mat(i,t-1),:)),1);
                            
                    
                            %Compute consumption 
                                if (a_min < a_mat(i,t)) && ( a_mat(i,t) < a_max)
                                    i_lo = min(n_a_fine-1, VFI_Toolbox.Grid_Inv(a_mat(i,t), n_a_fine, theta_a_f, a_min, a_max, "Poly" ));
                                    omega = min(1, max(0,(a_mat(i,t) - a_grid_fine(i_lo))/(a_grid_fine(i_lo+1) - a_grid_fine(i_lo)   )))  ;
                                    c_mat(i,t) = omega*G_c_fine(i_lo, eps_mat(i,t), age) + (1-omega)*G_c_fine(i_lo+1, eps_mat(i,t), age )  ;
                                elseif a_min==a_mat(i,t)
                                    c_mat(i,t) = G_c_fine(1,eps_mat(i,t), age);
                                else
                                    c_mat(i,t) = G_c_fine(end,eps_mat(i,t), age);
                                end
                                % Update age
                                h_mat(i,t)=1;
                        
                            end
                        end
                    end
                end
                
                %%%%%%%%%%%%%%%% Output %%%%%%%%%%%%%%%%
                

                PanelModel_Structure_Cohort.a_mat = a_mat ; % Output 1
                PanelModel_Structure_Cohort.c_mat = c_mat ; % Output 2
                PanelModel_Structure_Cohort.eps_mat = eps_mat ; % Output 3
                PanelModel_Structure_Cohort.h_mat = h_mat ; % Output 4

                Updated_PanelModel_Structure = PanelModel_Structure_Cohort;

            end % End of Simulate_Cohort function



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function Simulate_Panel_Timed

function ModelStructure_Output = Simulate_Panel_Timed(ModelStructure,PanelModel_Structure) % 
                
                % ModelStructure is the input structure (M in Julia code)
                % PanelModel_Structure is the structure of the panel parameters
                
                % ----------- Output variables ----------- 
                % Output1 = PM (updated version) 
                % Output1 is a Panel Model Structure 
                % ----------------------------------------  


                % Local Variables
                M = ModelStructure;
                PM = PanelModel_Structure;

                p = M.p;
                MP_eps = M.MP_eps;
                a_grid = M.a_grid;
                a_grid_fine = M.a_grid_fine;
                a_max = M.a_max;
                n_eps = M.n_eps;
                Gamma = M.Gamma;
                G_ap_fine = M.G_ap_fine;
                G_c_fine = M.G_c_fine;
                n_a_fine = M.n_a_fine;
                theta_a_f = M.theta_a_f;
                eps_grid = M.eps_grid;

                a_min = p.a_min;
                Surv_Pr = p.Surv_Pr;
                Max_Age = p.Max_Age;
                

                N_Panel = PM.N_Panel;
                T_Panel = PM.T_Panel;
                T_Simul = PM.T_Simul;
                t_mat =   PM.t_mat;
                
                % Timing starts! 
                time_now = tic();

                %Seed_Flag = 1;
                
                % Initialize Seed
                rng(PM.rng_seed);

                % Initialize cross-section vectors
                a_vec = zeros(N_Panel,1);
                c_vec = zeros(N_Panel,1);
                eps_vec = zeros(N_Panel,1);
                h_vec = zeros(N_Panel,1);

                % PDFs for epsilon and zeta
                Gamma_eps  = MP_eps.Pi;

                % Set median epsilon for newborns
                med_eps = int64(round(n_eps/2));

                % Set initial wealth to (close to) $1k
                b_ind = transpose(find(a_grid_fine>=1));
                b_ind = b_ind(1);
                b = a_grid_fine(b_ind);

                % Censor savings
                G_ap_fine = max(min(G_ap_fine,a_max), a_min);



                % Iterate Forward
                disp('Iterating Panel')

                for i=1:N_Panel
                    if mod(i,10000)==0
                        fprintf('   Dynasty %i \n', i)
                    end
                    % Draw initial conditions from stationary distribution (cross-section)
    
                    % Draw eps(i) and h(i) from (stationary) marginal
                    % distribution
                    eps_pdf = transpose(squeeze(sum(Gamma, [1,3])));
                    eps_vec = Functions_MonteCarlo.pdfrnd(eps_pdf, N_Panel);
                    h_pdf = squeeze(sum(M.Gamma, [1,2]));
                    h_vec = Functions_MonteCarlo.pdfrnd(h_pdf, N_Panel);
                    clear eps_pdf h_pdf

                    % Draw a(i) from conditional distribution
                    if h_vec(i)==1
                        a_vec(i)=b;
                    else
                        Gamma_aux = Gamma(:,eps_vec(i), h_vec(i))/sum(Gamma(:,eps_vec(i), h_vec(i)));
                        a_vec(i) = a_grid_fine(Functions_MonteCarlo.pdfrnd(Gamma_aux, 1));
                    end
    
                    for t=1:T_Simul
                        if mod(t,50)==0
                            fprintf(' Simulation Period %i \n', t)
                            fprintf(' E[a_t] = %4.2fk E[eps_t]=%3.3f E[h_t]=%4.1f, max(a_t) = %8.2f  \n', round(mean(a_vec),2), round(mean(eps_grid(eps_vec)),3), ...
                        round(mean(h_vec),1), round(max(a_vec),2) )
                        end
                        
                            Death_Draw = Functions_MonteCarlo.pdfrnd([Surv_Pr(h_vec(i)), 1-Surv_Pr(h_vec(i)) ],1);
    
                            if Death_Draw ==2 % Agent dies
                                a_vec(i)=b;
                                eps_vec(i)=med_eps;
                                h_vec(i)=1;
                            else % Agent lives
                                % Compute future assets: Linear interpolation (manually)
                                if (a_min<a_vec(i)) && (a_vec(i)<a_max)
                                    i_lo = min(n_a_fine-1, VFI_Toolbox.Grid_Inv(a_vec(i), n_a_fine, theta_a_f, a_min, a_max, "Poly" ));
                                    omega = min(1, max(0,(a_vec(i) - a_grid_fine(i_lo))/(a_grid_fine(i_lo+1) - a_grid_fine(i_lo))))  ;
                                    a_vec(i) = omega*G_ap_fine(i_lo, eps_vec(i), h_vec(i)) + (1-omega)*G_ap_fine(i_lo+1, eps_vec(i), h_vec(i))  ;
                                elseif a_vec(i)==a_max
                                     a_vec(i) = G_ap_fine(end, eps_vec(i), h_vec(i)) ;
                                elseif a_vec(i)==a_min
                                    a_vec(i) = G_ap_fine(1, eps_vec(i), h_vec(i)) ;
                                else
                                error('Error in simulation: assets not working')
                                end
                                a_vec(i) = min(a_max, max(a_min,a_vec(i)));
                                % Compute future eps(i)
                                eps_pdf = transpose(Gamma_eps(eps_vec(i),:));
                                eps_vec(i) = Functions_MonteCarlo.pdfrnd(eps_pdf,1);
                                clear eps_pdf
                                % Update age
                                h_vec(i) = h_vec(i)+1;
                            end
                            
                         
                            if t>= T_Simul - (T_Panel - 1)
                                %Compute consumption only for relevant periods 
                                if (a_min < a_vec(i)) && (a_vec(i) < a_max)
                                        i_lo = min(n_a_fine-1, VFI_Toolbox.Grid_Inv(a_vec(i), n_a_fine, theta_a_f, a_min, a_max, "Poly" ));
                                        omega = min(1, max(0,(a_vec(i) - a_grid_fine(i_lo))/(a_grid_fine(i_lo+1) - a_grid_fine(i_lo)   )))  ;
                                        c_vec(i) = omega*G_c_fine(i_lo, eps_vec(i), h_vec(i)) + (1-omega)*G_c_fine(i_lo+1, eps_vec(i), h_vec(i) )  ;
                                    elseif a_min==a_vec(i)
                                        c_vec(i) = G_c_fine(1,eps_vec(i), h_vec(i));
                                    else
                                        c_vec(i) = G_c_fine(end,eps_vec(i), h_vec(i));
                                end
                                 % Save results in panel
                                PM.a_mat(:, t-(T_Simul-T_Panel)) = a_vec ; % Assets
                                PM.eps_mat(:, t-(T_Simul-T_Panel)) = eps_vec ; % Assets
                                PM.h_mat(:, t-(T_Simul-T_Panel)) = h_vec ; % Assets
                                PM.c_mat(:, t-(T_Simul-T_Panel)) = c_vec ; % Assets
                            end
                    end
                    t_mat(i) = toc(time_now);
                    end
                PM.t_mat = t_mat;

                ModelStructure_Output = PM;
                return
               
end % End of Function Simulate_Panel_Timed



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function Simulate_Cohort_Timed


            function Updated_PanelModel_Structure = Simulate_Cohort_Timed(ModelStructure, PanelModel_Structure_Cohort, age_0) 
                
                
                % c_flag is a Boolean variable

                % ----------- Output contains 4 variables ----------- 
                % Output1 = a_mat
                % Output2 = c_mat
                % Output3 = eps_mat
                % Output4 = h_mat
                % Output5 = t_mat
                % ----------------------------------------  
                
                % Local Variables
                M = ModelStructure;
                M_P_cohort = PanelModel_Structure_Cohort;
                p = M.p;

                MP_eps = M.MP_eps;
                a_grid = M.a_grid;
                eps_grid = M.eps_grid;
                a_grid_fine = M.a_grid_fine;
                a_max = M.a_max;
                n_eps = M.n_eps;
                Gamma = M.Gamma;
                G_ap_fine = M.G_ap_fine;
                G_c_fine = M.G_c_fine;
                n_a_fine = M.n_a_fine;
                theta_a_f = M.theta_a_f;
                
                a_min = p.a_min;
                Surv_Pr = p.Surv_Pr;
                Max_Age = p.Max_Age;
                t_mat =   M_P_cohort.t_mat;


                N_Panel = M_P_cohort.N_Panel;
                

                % Max Cohort Age
                Max_C_Age = Max_Age - (age_0-1);

                % Initialize seed
                rng(M_P_cohort.rng_seed);

                % Initialize cross-section vector
                
                a_mat    = zeros([N_Panel, Max_C_Age]);
                c_mat    =  zeros([N_Panel, Max_C_Age]);
                eps_mat  = zeros([N_Panel, Max_C_Age]);
                h_mat    = zeros([N_Panel,Max_C_Age]);
                t_mat = zeros([N_Panel,1]);

                % PDF for epsilon
                Gamma_eps = MP_eps.Pi;

                % Set median epsilon for newborns
                med_eps = int64(round(n_eps/2));

                % Set initial wealth to (close to) $1k
                b_ind = transpose(find(a_grid_fine>=1));
                b_ind = b_ind(1);
                b = a_grid_fine(b_ind);

                % Censor savings
                G_ap_fine = max(min(G_ap_fine, a_max), a_min);
                

                % Timing starts! 
                time_now = tic();

                for i=1:N_Panel
                    if mod(i,10000)==0
                        fprintf('   Dynasty %i \n', i)
                    end
                    % Draw initial conditions from stationary distribution (cross-section)
    
                    % Draw eps(i) and h(i) from (stationary) marginal
                    % distribution
                    pdf_epsmat = transpose(sum(Gamma(:,:,age_0)/sum(sum(sum(Gamma(:,:,age_0)))),[1,3]));
                    eps_mat(:,1) = Functions_MonteCarlo.pdfrnd(pdf_epsmat,N_Panel);
                    h_mat(:,1)=1;
                    % Draw a(i) from conditional distribution
                    if age_0==1
                        a_mat(:,1) = b;
                        c_mat(:,1) = G_c_fine(b_ind,med_eps,1);
                    else
                        Gamma_aux = Gamma(:,eps_mat(i,1),age_0)/sum(Gamma(:,eps_mat(i,1), age_0));
                        a_ind = Functions_MonteCarlo.pdfrnd(Gamma_aux,1);
                        a_mat(i,1) = a_grid_fine(a_ind);
                        c_mat(i,1) = G_c_fine(a_ind, eps_mat(i,1), age_0);
                    end

                    for t=2:Max_C_Age
                        %if mod(t,5)==0
                        %ind = h_mat(:,t-1)>0; % Select only alive agents
                        %fprintf('Simulation Period %i \n',t)
                        %fprintf(' E[a_t] = %5.2fk E[eps_t]=%3.3f E[h_t]=%4.1f, max(a_t) = %8.2f  \n', round(mean(a_mat(ind,t-1)),2), round(mean(eps_grid(eps_mat(ind,t-1))),3), round(mean(h_mat(:,t-1)),1), round(max(a_mat(ind,t-1)),2))
                        %end
                        % Set current age of agents in cohort
                        age = t+(age_0-1);
                        if h_mat(i,t-1)==1
                            Death_Draw = Functions_MonteCarlo.pdfrnd([Surv_Pr(age), 1-Surv_Pr(age) ],1);
                            if Death_Draw ==2 % Agent dies, no need to update 
                                h_mat(i,t)=0;
                            else % Agent lives
                            % Compute future assets: Linear interpolation (manually)
                                if (a_min<a_mat(i,t-1)) && (a_mat(i,t-1)<a_max)
                                    i_lo = min(n_a_fine-1, VFI_Toolbox.Grid_Inv(a_mat(i,t-1), n_a_fine, theta_a_f, a_min, a_max, "Poly" ));
                                    omega = min(1, max(0,(a_mat(i,t-1) - a_grid_fine(i_lo))/(a_grid_fine(i_lo+1) - a_grid_fine(i_lo)   )))  ;
                                    a_mat(i,t) = omega*G_ap_fine(i_lo, eps_mat(i,t-1), age-1) + (1-omega)*G_ap_fine(i_lo+1, eps_mat(i,t-1), age-1 )  ;
                                elseif a_mat(i,t-1)==a_max
                                     a_mat(i,t) = G_ap_fine(end, eps_mat(i,t-1), age-1) ;
                                elseif a_mat(i,t-1)==a_min
                                    a_mat(i,t) = G_ap_fine(1, eps_mat(i,t-1),age-1) ;
                                else
                                error('Error in simulation: assets not working')
                                end
                            a_mat(i,t) = min(a_max, max(a_min,a_mat(i,t)));

                            % Compute future eps(i)
                            eps_mat(i,t) = Functions_MonteCarlo.pdfrnd(transpose(Gamma_eps(eps_mat(i,t-1),:)),1);
                            
                    
                            %Compute consumption 
                                if (a_min < a_mat(i,t)) && ( a_mat(i,t) < a_max)
                                    i_lo = min(n_a_fine-1, VFI_Toolbox.Grid_Inv(a_mat(i,t), n_a_fine, theta_a_f, a_min, a_max, "Poly" ));
                                    omega = min(1, max(0,(a_mat(i,t) - a_grid_fine(i_lo))/(a_grid_fine(i_lo+1) - a_grid_fine(i_lo)   )))  ;
                                    c_mat(i,t) = omega*G_c_fine(i_lo, eps_mat(i,t), age) + (1-omega)*G_c_fine(i_lo+1, eps_mat(i,t), age )  ;
                                elseif a_min==a_mat(i,t)
                                    c_mat(i,t) = G_c_fine(1,eps_mat(i,t), age);
                                else
                                    c_mat(i,t) = G_c_fine(end,eps_mat(i,t), age);
                                end
                                % Update age
                                h_mat(i,t)=1;
                        
                            end
                        end


                    end % Age(t)
                    t_mat(i) = toc(time_now); % End of loop for i
                end
                PanelModel_Structure_Cohort.a_mat = a_mat ; % Output 1
                PanelModel_Structure_Cohort.c_mat = c_mat ; % Output 2
                PanelModel_Structure_Cohort.eps_mat = eps_mat ; % Output 3
                PanelModel_Structure_Cohort.h_mat = h_mat ; % Output 4
                PanelModel_Structure_Cohort.t_mat = t_mat ; % Output 5

                
                %%%%%%%%%%%%%%%% Output %%%%%%%%%%%%%%%%                

                Updated_PanelModel_Structure = PanelModel_Structure_Cohort;

            end % End of Simulate_Cohort_Timed function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end % End of methods(static)

    end % end of Classdef