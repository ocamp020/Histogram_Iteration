
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
                M_P_cohort = PanelModel_Structure;
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
                b = transpose(a_grid(ind_inner));

                % Censor savings
                G_ap_fine = max(min(G_ap_fine, a_max), a_min);
                
                % Draw initial conditions from stationary distribution (cross-section)
                disp("Initializing Simulation")

                % Draw epsilon(i) from (stationary) marignal distribution
                pdf_epsmat = sum(Gamma(:,:,age_0)/sum(sum(sum(Gamma(:,:,age_0)))),[1,3])
                eps_mat(:,1) = pdfrnd(pdf_epsmat,N_Panel)
                % ------------------ Line 213 from Julia! 

                %%%%%%%%%%%%%%%% Output %%%%%%%%%%%%%%%%
                

                PanelModel_Structure_Cohort.a_mat = a_mat ; % Output 1
                PanelModel_Structure_Cohort.c_mat = c_mat ; % Output 2
                PanelModel_Structure_Cohort.eps_mat = eps_mat ; % Output 3
                PanelModel_Structure_Cohort.h_mat = h_mat ; % Output 4

                Updated_PanelModel_Structure = PanelModel_Structure_Cohort;

            end % End of Simulate_Cohort function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Panel Simulation
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


                %M = M_Aiyagari;
                %PM = M_P;

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

                % Set initial wealth to(close to) $1k

                % Set initial wealth to (close to) $1k
                b_ind = transpose(find(a_grid_fine>=1));
                b = transpose(a_grid(ind_inner));

                % Censor savings
                % G_ap_fine = max(min(G_ap_fine,a_max), a_min);

                % Draw initial conditions from stationary distribution (cross-section)
                fprintf(' Initializing simulation \n')

                % Draw eps(i) and zeta(i) from (stationary) marginal
                % distribution
                eps_vec = Functions_MonteCarlo.pdfrnd(MP_eps.PDF, N_Panel);
                zeta_vec = Functions_MonteCarlo.pdfrnd(MP_zeta.PDF, N_Panel);

                if Seed_Flag==1 % Use conditional distribution in Γ from histogram method 
                    for i=1:N_Panel
                        Gamma_aux = Gamma(:,eps_vec(i), zeta_vec(i))/sum(Gamma(:,eps_vec(i), zeta_vec(i)));
                        a_vec(i) = a_grid_fine(Functions_MonteCarlo.pdfrnd(Gamma_aux, 1));
                    end
                else
                    a_vec = unifrnd(a_min,a_grid_fine(n_cut_fine),N_Panel);
                end

                
                %*************************************************************
                %*************************************************************
                % Importing julia's a_vec, eps_vec, zeta_vec 
                % bc they are random vectors so I just want to check the
                % code has the same output given that it 
                % works with the same input numbers
                %a_vec = readmatrix('/Users/cyberdim/Dropbox/WESTERN_ECONOMICS/RA_Baxter/Matlab_Hist_Iter/Infinitely_Lived/File_Folder/a_vec_Julia.csv');
                %eps_vec = readmatrix('/Users/cyberdim/Dropbox/WESTERN_ECONOMICS/RA_Baxter/Matlab_Hist_Iter/Infinitely_Lived/File_Folder/eps_vec_Julia.csv');
                %zeta_vec = readmatrix('/Users/cyberdim/Dropbox/WESTERN_ECONOMICS/RA_Baxter/Matlab_Hist_Iter/Infinitely_Lived/File_Folder/zeta_vec_Julia.csv');
                %*************************************************************
                %*************************************************************



                fprintf('E[a_0] = $%8.2f k , E[eps_0]= %8.3f , E[zeta_0] = %8.3f, max(a_0) = %8.2f \n', ...
                    round(mean(a_vec),2), round(mean(eps_grid(eps_vec)),3), ...
                    round(mean(zeta_grid(zeta_vec)),3), round(max(a_vec),2)  )
                
                % Initialize moments
                pct_list = [50; 90;99;99.9;99.99];
                pct_vec = prctile(a_vec, pct_list)   ; % This percentiles even for the same a_vec differ between Julia and Matlab!
                
                ts_vec = zeros(5,1);
                for p =1:numel(pct_list) 
                    pct_vec_inner = pct_vec(p);
                    I = find(a_vec>=pct_vec_inner);
                    ts_vec(p) = 100*sum(a_vec(I))/sum(a_vec);
                end
                moment_vec_new = horzcat(mean(a_vec), transpose(pct_vec(:)));
                moment_vec_old = moment_vec_new;
                Simul_dist=10;

                % Iterate forward
                disp("Iterating Panel")
                t=0;
                while (t<=T_Simul)||(Simul_dist>Simul_tol)
                    t = t+1;
                %for t=1:T_Simul
                    if mod(t,50)==0
                        fprintf('Simulation period %i \n', t)
                        fprintf('E[a_t] = $%8.2f k , E[eps_t]= %5.3f , E[zeta_t] = %5.3f, max(a_t) = %8.4f , dist = %8.4f \n', ...
                    round(mean(a_vec),2), round(mean(eps_grid(eps_vec)),3), ...
                    round(mean(zeta_grid(zeta_vec)),3), round(max(a_vec),2), round(Simul_dist,4)  )
                    end


                    % Simulate each dinasty
                    for i=1:N_Panel
                        % Compute future assets: Linear interpolation (manually)
                        if (a_min<a_vec(i)) && (a_vec(i)<a_max)
                            i_lo = min(n_a_fine-1, VFI_Toolbox.Grid_Inv(a_vec(i), n_a_fine, theta_a_f, a_min, a_max, "Poly" ));
                            omega = 1 - min(1, max(0,(a_vec(i) - a_grid_fine(i_lo))/(a_grid_fine(i_lo+1) - a_grid_fine(i_lo)   )))  ;
                            a_vec(i) = omega*G_ap_fine(i_lo, eps_vec(i), zeta_vec(i)) + (1-omega)*G_ap_fine(i_lo+1, eps_vec(i), zeta_vec(i) )  ;
                        elseif a_vec(i)==a_max
                            a_vec(i) = G_ap(end, eps_vec(i), zeta_vec(i)) ;
                        elseif a_vec(i)==a_min
                            a_vec(i) = G_ap(1, eps_vec(i), zeta_vec(i)) ;
                        else
                            error('Error in simulation: assets not working')
                        end
                        % Censor assets to be between a_min and a_max
                        a_vec(i) = min(a_max, max(a_min, a_vec(i)))   ;

                        % Compute future ϵ[i] and ζ[i]
                        eps_vec(i)  = Functions_MonteCarlo.pdfrnd(Gamma_eps (eps_vec(i),:),1) ;
                        zeta_vec(i) = Functions_MonteCarlo.pdfrnd(Gamma_zeta(zeta_vec(i),:),1);

                        % Compute consumption only for relevant periods 
                        %if t>=T_Simul - (T_Panel - 1 )
                            %if (a_min < a_vec(i)) && (a_vec(i) < a_max)
                                %G_c_ip = griddedInterpolant(a_grid, G_c(:, eps_vec(i), zeta_vec(i)), 'spline', 'nearest');
                                %c_vec(i) = G_c_ip(a_vec(i));
                            %elseif a_min == a_vec(i)
                            %    c_vec(i) = G_c(1, eps_vec(i), zeta_vec(i));
                            %else
                            %    c_vec(i) = G_c(end, eps_vec(i), zeta_vec(i));
                            %end
                        %end
                    end

                        % Save results in panel
                        if t >= (T_Simul - (T_Panel - 1))
                            a_panel = horzcat(a_panel(:,2:end), a_vec);
                            eps_panel = horzcat(eps_panel(:,2:end), eps_vec);
                            zeta_panel = horzcat(zeta_panel(:,2:end), zeta_vec) ;
                        end

                        %if t>+ T_Simul - (T_Panel - 1)
                        %    M_P.a_mat(:, t-(T_Simul-T_Panel)) = a_vec ; % Assets
                        %    M_P.eps_mat(:, t-(T_Simul-T_Panel)) = eps_vec ; % Assets
                        %    M_P.zeta_mat(:, t-(T_Simul-T_Panel)) = zeta_vec ; % Assets
                        %    M_P.c_mat(:, t-(T_Simul-T_Panel)) = c_vec ; % Assets
                        %end

                        % Compute moments
                        pct_list = [50; 90;99;99.9;99.99];
                        pct_vec = prctile(a_vec, pct_list);
                        ts_vec = zeros(5,1);
                        for p=1:numel(pct_list)
                            pct_vec_inner = pct_vec(p);
                            I = find(a_vec>=pct_vec_inner);
                            ts_vec(p) = 100*sum(a_vec(I))/sum(a_vec);
                        end
                        moment_vec_new = horzcat(mean(a_vec), transpose(pct_vec(:)) );
                        % Compute moments distance
                        Simul_dist = max(max(max(abs((moment_vec_new./moment_vec_old)-1))));
                        moment_vec_old = moment_vec_new;


                end
                fprintf(" Iteration completed after %i periods with dist=%8.4f \n",t, Simul_dist )
                % Save results in panel
                PM.a_mat = a_panel; % Assets
                PM.eps_mat = eps_panel ; % Labor efficiency
                PM.zeta_mat = zeta_panel; % Returns

                % Compute Consumption
                fprintf(" Computing consumption for panel of %i periods  \n ", T_Panel)
                for t=1:T_Panel

                    %Extract States
                    a_vec = a_panel(:,t);
                    eps_vec = eps_panel(:,t);
                    zeta_vec = zeta_panel(:,t);

                    % Compute consumption for current period
                    for i=1:N_Panel
                        if (a_min < a_vec(i)) && (a_vec(i) < a_max)
                            x = a_grid;
                            y = G_c(:,eps_vec(i),zeta_vec(i));
                            sp = fn2fm(spline(x,y),'B-');
                            yy = fnval(sp, a_vec(i));
                            c_vec(i) = transpose(yy(1,1)); 
                        elseif a_min==a_vec(i)
                            c_vec(i) = G_c(1,eps_vec(i),zeta_vec(i));
                        else
                            c_vec(i) = G_c(end, eps_vec(i), zeta_vec(i))  ;
                        end
                    end

                    PM.c_mat(:,t) = c_vec; % Consumption

                end

                ModelStructure_Output = PM;
                return
               
            end % End of Function Simulate_Panel




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Simulate the panel - Dinasty by Dinasty


            function ModelStructure_Output = Simulate_Panel_Dynasty(ModelStructure, PanelModel_Structure)
                    
                % ModelStructure is the input structure
                % PanelModel_Structure is the structure of the panel parameters
           
                
                
                % Local Variables
                M = ModelStructurel;
                PM = PanelModel_Structure;
                p = M.p;
                MP_eps = M.MP_eps;
                MP_zeta = M.MP_zeta;
                a_grid = M.a_grid;
                a_grid_fine = M.a_grid_fine;
                a_max = M.a_max;
                Gamma = M.Gamma;
                G_ap = M.G_ap;
                G_ap_fine = M.G_ap_fine;
                G_c = M.G_c;
                n_a_fine = M.n_a_fine;
                theta_a_f = M.theta_a_f;

                a_min = M.p.a_min;

                N_Panel = PM.N_Panel;
                T_Panel = PM.T_Panel;
                T_Simul = PM.T_Simul;
                
                % PDFs for epsilon and zeta
                Gamma_eps = M.MP_eps.Pi;
                Gamma_zeta = M.MP_zeta.Pi;

                % Initialize Seed
                rng(PM.rng_seed);

                % Censor Savings
                G_ap = max(min(G_ap,a_max),a_min);
                G_ap_fine = max(min(G_ap_fine,a_max),a_min);

                % Iterate forward 
                disp(' Iterating forward')
                for i=1:N_Panel
                    tic 
                    
                    if mod(i,1000)==0
                        fprintf(' Dynasty %i \n', i)
                    end

                    % Draw initial conditions from stationary distribution (cross-section)
                    %  Draw ϵ[i] and ζ[i] from (stationary) marginal distribution 

                    eps_vec = Functions_MonteCarlo.pdfrnd(MP_eps.PDF,1);
                    zeta_vec = Functions_MonteCarlo.pdfrnd(MP_zeta.PDF,1);
                    % Draw a[i] from conditional distribution 
                    Gamma_aux = Gamma(:, eps_vec, zeta_vec)/sum(Gamma(:,eps_vec, zeta_vec));
                    a_vec = a_grid_fine(Functions_MonteCarlo.pdfrnd(Gamma_aux,1));

                    for t=1:T_Simul

                        % if mod(t,100)==0
                        % fprintf(' Simulation Period %i \n', t)
                        % end

                        % Compute future assets: Linear interpolation
                        % (manually)

                        if (a_min<a_vec) && (a_vec<a_max)
                            i_lo = min(n_a_fine-1, VFI_Toolbox.Grid_Inv(a_vec, n_a_fine, theta_a_f, a_min, a_max, "Poly"  ));
                            omega = min(1, max(0, (a_vec - a_grid_fine(i_lo))/(a_grid_fine(i_lo+1)-a_grid_fine(i_lo) ) ));
                            a_vec = omega*G_ap_fine(i_lo,eps_vec,zeta_vec) + (1-omega)*G_ap_fine(i_lo+1, eps_vec, zeta_vec);
                        elseif a_vec == a_max 
                            a_vec = G_ap(end, eps_vec, zeta_vec);
                        elseif a_vec == a_min
                            a_vec =  G_ap(1, eps_vec, zeta_vec);
                        else
                            disp('Error in simulation: assets not working')
                        end
                        a_vec = min(a_max, max(a_min, a_vec));

                        % Compute future eps(i) and zeta(i)

                        eps_vec = Functions_MonteCarlo.pdfrnd(Gamma_eps(eps_vec,:),1);
                        zeta_vec = Functions_MonteCarlo.pdfrnd(Gamma_zeta(zeta_vec,:),1);

                        % Compute consumption only for relevant periods

                        if t>T_Simul-(T_Panel-1)
                            if (a_min < a_vec) && (a_vec < a_max)
                                x = a_grid;
                                y = G_c(:,eps_vec,zeta_vec);
                                sp = fn2fm(spline(x,y),'B-');
                                yy = fnval(sp, a_vec);
                                c_vec = transpose(yy); 
                            elseif a_min==a_vec
                                c_vec = G_c(1,eps_vec,zeta_vec);
                            else
                                c_vec = G_c(end, eps_vec, zeta_vec);
                            end
                        end

                        % Save results in panel
                        if t>T_Simul-(T_Panel-1)
                            PM.a_mat(i,t-(T_Simul-T_Panel))=a_vec; % Assets
                            PM.eps_mat(i,t-(T_Simul-T_Panel))=eps_vec; % Labor efficiency
                            PM.zeta_mat(i,t-(T_Simul-T_Panel))=zeta_vec; % Returns
                            PM.c_mat(i,t-(T_Simul-T_Panel))=c_vec; %  Consumption
                        end
                    end

                    toc
                end 


                ModelStructure_Output = PM;
                return

            end % End of Function Simulate_Panel_Dynasty
 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




        end % End of methods(static)

    end % end of Classdef