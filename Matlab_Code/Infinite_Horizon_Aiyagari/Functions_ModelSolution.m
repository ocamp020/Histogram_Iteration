% Functions
classdef Functions_ModelSolution
    methods(Static)

    % Utility function
        function util = utility(c,par)
            if par.sigma > 1
                util = (c.^(1-par.sigma))/(1-par.sigma);
            else
                util = log(c);
            end
        end

    % Derivative of the utility function
        function derivative = d_utility(c,par)
            derivative = c.^(-par.sigma);
        end

    % Inverse of the derivative
        function duinv = d_utility_inv(x,par)
            duinv = x.^(-1/par.sigma);
        end
        
% ------------------------------------------------------------------------------------------- 
% ------------------------------------------------------------------------------------------- 
% -------------------------------------------------------------------------------------------

    % Parameter structure
    function par_structure = Parameters(beta,sigma, ...
            rho_eps,sigma_eps,rho_zeta,sigma_zeta, ...
            r,w,a_min,max_iter,dist_tol,dist_tol_Delta, ...
            eta, Hist_max_iter, Hist_tol,Hist_eta, c_min)

            par_structure = struct('beta',beta, 'sigma',sigma, 'rho_eps', rho_eps, ...
                'sigma_eps',sigma_eps,'rho_zeta',rho_zeta, 'sigma_zeta',sigma_zeta,...
    'r',r, 'w',w,...
    'a_min',a_min,...
    'max_iter',max_iter,...
    'dist_tol',dist_tol,'dist_tol_Delta',dist_tol_Delta,'eta',eta,...
    'Hist_max_iter',Hist_max_iter,'Hist_tol',Hist_tol,'Hist_eta',Hist_eta,...
    'c_min',c_min); 
            return
    end % End of function Parameters


    function updated_par_structure = par2modify(parameter_structure, parameter_requested, parameter_new_value)
        % This function checks if the requested parameter is in the
        % parameter structure. If so, requests the function 'Parameters' for 
        % an updated parameter structure. If the requested parameter is not
        % there, displays an error sign.

        % This is what I though to use to overcome the @kw macro/structure
        % in Julia

        % The requested parameter should be written as: 'parameter_requested' .
        
        old = parameter_structure;
        if any(strcmp(fieldnames(parameter_structure),parameter_requested))
            old.(sprintf(parameter_requested)) = parameter_new_value;
            values = struct2cell(old);
            %fields_inner
            updated_par_structure = Functions_ModelSolution.Parameters(values{:});
        else
            error('Parameter not in the structure')
        end
        clear old
        return
        
    end % End of function par2modify

        
% ------------------------------------------------------------------------------------------- 
% -------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------


    % Model structure

    function model_structure = ModelStructure(par_structure, a_max, theta_a, theta_a_f, n_a ,  n_a_fine, n_zeta, n_eps, method, read_flag)
        % Inputs: parameter structure plus 9 scalar inputs that allow the function to create the
        % model structure
        
        % method: 1 for Kronecker and 2 for loops in expectation of PFI 
            % read_flag: Boolean for reading results from file. 
            model_structure = struct('p', par_structure, ...
                'a_max', a_max, 'theta_a', theta_a, ...
                'theta_a_f', theta_a_f, 'n_a',n_a, 'n_a_fine',n_a_fine, ...
                'n_zeta',n_zeta, 'n_eps', n_eps, ...
                'method', method , 'read_flag', read_flag);
            par = model_structure.p;
            % model_structure complement
                model_structure.a_grid = VFI_Toolbox.Make_Grid(model_structure.n_a, model_structure.theta_a, par.a_min, model_structure.a_max,"Poly");
                model_structure.a_grid_fine = VFI_Toolbox.Make_Grid(model_structure.n_a_fine, model_structure.theta_a_f, par.a_min, model_structure.a_max,"Poly");
                model_structure.n_cut_fine = double(VFI_Toolbox.Grid_Inv(20000,model_structure.n_a_fine,model_structure.theta_a_f, par.a_min, model_structure.a_max,"Poly"));
                % Interest rate process
                model_structure.MP_zeta = struct('N',{}, 'grid',{}, 'Pi',{}, 'PDF',{},'CDF',{});
                [model_structure.MP_zeta(1).N,model_structure.MP_zeta(1).grid,model_structure.MP_zeta(1).Pi,model_structure.MP_zeta(1).PDF,model_structure.MP_zeta(1).CDF] = VFI_Toolbox.Tauchen86(par.rho_zeta,par.sigma_zeta,model_structure.n_zeta,1.96);
                model_structure.zeta_ref = 1/dot(exp(model_structure.MP_zeta(1).grid),model_structure.MP_zeta(1).PDF);
                model_structure.zeta_grid = model_structure.zeta_ref.*exp(model_structure.MP_zeta(1).grid);
                % Labor productivity process
                model_structure.MP_eps = struct('N',{}, 'grid',{}, 'Pi',{}, 'PDF',{},'CDF',{});
                [model_structure.MP_eps(1).N,model_structure.MP_eps(1).grid,model_structure.MP_eps(1).Pi,model_structure.MP_eps(1).PDF,model_structure.MP_eps(1).CDF] = VFI_Toolbox.Rouwenhorst95(par.rho_eps,par.sigma_eps,model_structure.n_eps);
                model_structure.eps_ref = 1/dot(exp(model_structure.MP_eps(1).grid),model_structure.MP_eps(1).PDF);
                model_structure.eps_grid = model_structure.eps_ref.*exp(model_structure.MP_eps(1).grid);
                % State matrices
                model_structure.a_mat = repmat(transpose(model_structure.a_grid),[1,model_structure.n_eps,model_structure.n_zeta]);
                model_structure.eps_mat = repmat(model_structure.eps_grid,[model_structure.n_a,1,model_structure.n_zeta]);
                model_structure.zeta_mat = repmat(reshape(model_structure.zeta_grid,[1,1,model_structure.n_zeta]),model_structure.n_a, model_structure.n_eps,1);
                model_structure.a_mat_fine =  repmat(transpose(model_structure.a_grid_fine),[1,model_structure.n_eps,model_structure.n_zeta]);
                model_structure.eps_mat_fine = repmat(model_structure.eps_grid,[model_structure.n_a_fine,1,model_structure.n_zeta]);
                model_structure.zeta_mat_fine = repmat(reshape(model_structure.zeta_grid,[1,1,model_structure.n_zeta]),model_structure.n_a_fine, model_structure.n_eps,1);
                % Value and policy functions 
                model_structure.V = zeros([model_structure.n_a, model_structure.n_eps, model_structure.n_zeta]); % Value function
                model_structure.G_ap = zeros([model_structure.n_a, model_structure.n_eps, model_structure.n_zeta]); % Policy Function for Capital/Assets
                model_structure.G_c = zeros([model_structure.n_a, model_structure.n_eps, model_structure.n_zeta]); % Policy Function for Consumption
                model_structure.V_fine = zeros([model_structure.n_a_fine, model_structure.n_eps, model_structure.n_zeta]); % Value function
                model_structure.G_ap_fine = zeros([model_structure.n_a_fine, model_structure.n_eps, model_structure.n_zeta]); % Policy Function for Capital/Assets
                model_structure.G_c_fine = zeros([model_structure.n_a_fine, model_structure.n_eps, model_structure.n_zeta]); % Policy Function for Consumption
                % Distribution
                model_structure.Gamma = (1/(model_structure.n_cut_fine*model_structure.n_eps*model_structure.n_zeta)).*cat(1,ones(model_structure.n_cut_fine,model_structure.n_eps,model_structure.n_zeta), zeros(model_structure.n_a_fine-model_structure.n_cut_fine, model_structure.n_eps, model_structure.n_zeta));
                model_structure.H_ind = rand(model_structure.n_a_fine, model_structure.n_eps, model_structure.n_zeta); % Index for discretization of savings choice  
                model_structure.H_omega_lo = rand(model_structure.n_a_fine, model_structure.n_eps, model_structure.n_zeta, model_structure.n_eps, model_structure.n_zeta); % Index for discretization of savings choice  
                model_structure.H_omega_hi = rand(model_structure.n_a_fine, model_structure.n_eps, model_structure.n_zeta, model_structure.n_eps, model_structure.n_zeta);

            return


    end % End of ModelStructure

      % ***************************************************
   
 function updated_model_structure = model2modify(model_structure, parameter_requested, parameter_new_value)
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
            updated_model_structure = Functions_ModelSolution.ModelStructure(values{1:10});
        else
            error('Parameter not in the structure')
        end
        clear old
        return
        
    end % End of function model2modify


        
% ------------------------------------------------------------------------------------------- 
% -------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------
            
        % Policy Function Iteration
        function ModelStructure_Output = PFI_Fixed_Point(FT,ModelStructure, PolicyFunctionMatrix) % 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % When calling this function it should be like:
            % ModelStructure_Output = PFI_Fixed_Point(FT,ModelStructure,PolicyFunctionMatrix)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % FT is the function type (i.e. the type of Bellman operator)
            % ModelStructure is the structure keeping the model values, including the parameters 
            % PolicyFunctionMatrix: is the matrix for assets (G_ap_old in the Julia code)

            M = ModelStructure; % To avoid typing everywhere ModelStructure
            

            % The way the Julia code is written, PolicyFunctionMatrix = G_ap_old
            % is forced to be empty. So to make it work equivalently before
            % raising this with Baxter and Sergio, I make it an empty
            % function here always (even when read as a filled matrix).

            PolicyFunctionMatrix = zeros(0,0,1);
            % Initialize variables for loop
            if isempty(PolicyFunctionMatrix) == 1
                G_ap_old = (M.p.r*M.zeta_mat+1).*M.a_mat;
            end

            G_dist_new = 100    ; % Initialize distance
            G_dist_old = 1      ; % Initialize distance
            G_dist_change = 1   ; % Initialize distance
            fprintf('-------------------------------- \n');
            fprintf('-------------------------------- \n');
            fprintf('PFI - n_eps = %d, n_zeta=%d, n_a = %d, theta_a =%d, r=%d \n', M.n_eps, M.n_zeta, M.n_a, M.theta_a, M.p.r);
            for iter= 1:M.p.max_iter
                % Update old distance and iterations
                G_dist_old = G_dist_new;
                % Update value function
                if FT=='T_EGM_G'
                    M.G_ap=G_ap_old;
                    [G_ap_new, G_c] = Functions_ModelSolution.T_EGM_G(M);
                end
                % Update new distance and iterations
                G_dist_new = sqrt(norm(G_ap_new-G_ap_old,"fro"));
                %fprintf('----------------------------------- \n')
                %fprintf('G_dist_new %d \n',G_dist_new )
                %fprintf('----------------------------------- \n')
                % Update change in distance
                G_dist_change = abs(G_dist_new-G_dist_old);
                %fprintf('----------------------------------- \n')
                %fprintf('G_dist_change %d \n',G_dist_change )
                %fprintf('----------------------------------- \n')
                % Update old function
                G_ap_old = (1-M.p.eta)*G_ap_new + M.p.eta*G_ap_old;
                
                %fprintf('-------------------------------- \n');
                %fprintf('PFI Loop: iter = %d, dist = %.4d \n',iter,G_dist_new)
                %fprintf('-------------------------------- \n');
                
                % Update change in convergence criteria
                % Report progress
                if  mod(iter,250)==0
                    fprintf('PFI Loop: iter = %d, dist = %.4d \n',iter,G_dist_new)
                end
                
                % Check convergence and return results

                if G_dist_new <= M.p.dist_tol || ((G_dist_change<=M.p.dist_tol_Delta)&&(G_dist_new<=M.p.dist_tol*10))
                    if G_dist_new<=M.p.dist_tol
                        fprintf("Distance converged: Iterations = %d, Distance = %.4d \n",iter,G_dist_new)
                    elseif ((G_dist_change<=M.p.dist_tol_Delta)&&(G_dist_new<=M.p.dist_tol*10))
                        fprintf("Distance fluctuating: Iterations = %d, Distance = %.4d \n", iter, G_dist_new)
                        fprintf("Change in distance converged: Iterations = %d, Δ Distance = %.4d \n",iter, G_dist_change)
                    end
                    disp("---------------------------------")
                    % Check borrowing constraint
                    if any(G_ap_new<M.p.a_min)
                        error('Borrowing Constraint Violated')
                    end
                    % Interpolate to fine grid
                    G_ap_fine = zeros(M.n_a_fine,M.n_eps,M.n_zeta);
                    G_c_fine = zeros(M.n_a_fine,M.n_eps,M.n_zeta);
                    for iter_eps=1:M.n_eps
                    for iter_zeta=1:M.n_zeta
                            % ------------------------- ------------------------- -------------------------
                            % For this part I had to do a slight
                            % modification to employ the spline function in
                            % Matlab. 
                            x = M.a_grid;
                            y = G_ap_new(:,iter_eps,iter_zeta);
                            sp = fn2fm(spline(x,y),'B-');
                            yy = fnval(sp, M.a_grid_fine);
                            G_ap_fine(:,iter_eps,iter_zeta) = transpose(yy); % Trying to get a B-spline in Matlab this way does not seem to solve the numerical differences.

                            %A_ip =  griddedInterpolant(M.a_grid, G_ap_new(:,iter_eps,iter_zeta), 'spline', 'previous');
                            %G_ap_fine(:,iter_eps,iter_zeta) = A_ip(M.a_grid_fine);

                            clear y yy sp

                            % Now consumption

                            y = G_c(:,iter_eps, iter_zeta);
                            sp = fn2fm(spline(x,y), 'B-');
                            yy = fnval(sp, M.a_grid_fine);
                            G_c_fine(:,iter_eps,iter_zeta) = transpose(yy);
                            %C_ip =  griddedInterpolant(M.a_grid, G_c(:,iter_eps, iter_zeta), 'spline', 'previous');
                            %G_c_fine(:,iter_eps,iter_zeta) = C_ip(M.a_grid_fine);

                            clear y yy sp
                            
                    end
                    end
                    % Update model
                    M.G_ap = G_ap_new;
                    M.G_c = G_c;
                    M.G_ap_fine = G_ap_fine;
                    M.G_c_fine = G_c_fine ;
                    ModelStructure_Output = M;
                    return
                end
            end
            error('Error in PFI - Solution not found')
        end % End of Function
        

% ------------------------------------------------------------------------------------------- 
% -------------------------------------------------------------------------------------------
        
        % Bellman operator - EGM - Iterate on Policy Functions


        function [Output1, Output2] = T_EGM_G(Model_Structure) 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % When calling this function it should be like:
            % [M.G_ap, M.G_c] = TT_EGM_G(M) 
            % Where M is the name of the Model Structure.
            % Output1 is G_ap
            % Output2 is G_c
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % p, n_ζ, n_ϵ, MP_ϵ, MP_ζ, n_a, G_ap, method = M
            % β, a_min, r, w = p
            
            M = Model_Structure;
            % Local variables
            p = M.p;
            n_zeta = M.n_zeta;
            n_eps = M.n_eps;
            MP_eps = M.MP_eps;
            MP_zeta = M.MP_zeta;
            n_a = M.n_a;
            G_ap = M.G_ap;
            method = M.method;
            zeta_mat = M.zeta_mat;
            a_mat= M.a_mat;
            eps_mat = M.eps_mat;
            a_grid=M.a_grid;
            eps_grid = M.eps_grid;

            beta = M.p.beta;
            a_min = M.p.a_min;
            r = M.p.r;
            w = M.p.w;
            
            % Define RHS of Euler equation for each (a',ϵ,ζ)
            % First dim is tomorrow's a in a fixed grid, second dim is present ϵ, third dim is present ζ
            
            if method==1
                Pi_joint = kron(MP_zeta.Pi, MP_eps.Pi);
                Euler_RHS = beta*Pi_joint*(reshape((r*zeta_mat+1).* ...
                    Functions_ModelSolution.d_utility((r*zeta_mat+1).*a_mat + w*eps_mat - G_ap, p ), [n_a, n_eps*n_zeta])');
                Euler_RHS = reshape(Euler_RHS', [n_a, n_eps, n_zeta]);
            elseif method==2
                % Loop Calculation
                aux_1 = (r*zeta_mat + 1 ).*Functions_ModelSolution.d_utility((r*zeta_mat+1).*a_mat + w*eps_mat - G_ap  ,p);
                aux_2 = randi([0,1],[n_eps,n_zeta]);
                Euler_RHS = randi([0,1],[n_a,n_eps,n_zeta]);
                for i_zeta=1:n_zeta % Today's zeta
                    Pr_zetap = MP_zeta.Pi(i_zeta,:);
                    for i_eps = 1:n_eps % Today's epsilon
                        Pr_epsp = MP_eps.Pi(i_eps,:);
                        for i_ap=1:n_a % Tomorrow's assets
                            for i_zetap=1:n_zeta
                                for i_epsp=1:n_eps
                                    aux_2(i_epsp,i_zetap) = beta*aux_1(i_ap,i_epsp,i_zetap)*Pr_epsp(i_epsp)*Pr_zetap(i_zetap);
                                end
                            end
                            Euler_RHS(i_ap,i_eps,i_zeta) = sum(aux_2(:,:));
                        end
                    end
                end
            end
            % Check monotonicity
            if any(Euler_RHS<0)
                error("RHS must be monotone for EGM to work")
            end
            % Define consumption from Euler equation
            C_endo = Functions_ModelSolution.d_utility_inv(Euler_RHS,p);
            % Define endogenous grid on assets
            A_endo = (C_endo + a_mat - w*eps_mat)./(r*zeta_mat+1);
            % Interpolate functions on exogenous grid
            G_c = zeros([n_a,n_eps,n_zeta]);
            for i_zeta=1:n_zeta
                for i_eps=1:n_eps
                    % Sort A_endo for interpolation
                    [A_aux,sort_ind] = sort(A_endo(:,i_eps,i_zeta));
                    C_aux = C_endo(:,i_eps,i_zeta);
                    C_aux = C_aux(sort_ind);
                    % Check boundary condition
                    if min(A_aux)>a_min
                        a_vec = transpose(a_grid(a_grid<min(A_aux)));
                        A_aux = [a_vec ; A_aux];
                        C_aux = [ ((r*zeta_mat(i_zeta)+1).*a_vec + w*eps_grid(i_eps)  - a_min) ; C_aux];
                    end
                    %C_ip = griddedInterpolant(A_aux, C_aux,'spline') ;   % 
                    C_ip = griddedInterpolant(A_aux, C_aux,'spline','nearest') ; % This is the closest to Julia Dierckx.jl interpolation routine.
                    G_c(:,i_eps,i_zeta) = C_ip(a_grid);
                    Ap_aux = (r*zeta_mat(i_zeta)+1).*transpose(a_grid) + w*eps_grid(i_eps) - G_c(:,i_eps,i_zeta);
                end
            end
            % Update policy function
            G_ap = (r*zeta_mat + 1).*a_mat + w*eps_mat - G_c;
            % Adjust for numerical error
            ind = find(abs(G_ap - a_min)<= 1e-10 );
            G_ap(ind) = a_min;
            G_c(ind) = (r*zeta_mat(i_zeta)+1)*a_mat(ind) + w*eps_mat(ind) - a_min;
            % Check for borrowing constraint
            % if any(G_ap<a_min)
                %error('Borrowing Constraint Violated')
            % end 
            % Return Results
            Output1 = G_ap;
            Output2 = G_c;
        end % End of T_EGM_G function

 
% ------------------------------------------------------------------------------------------- 
% -------------------------------------------------------------------------------------------
 
% Histogram method loop
% This function updates 

        function ModelStructure_Output = Histogram_Method_Loop(ModelStructure, N_H, Gamma_0)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % When calling this function it should be like:
            % ModelStructure_Output = Histogram_Method_Loop(M, N_H, Gamma_0) 
            % The outputs of this function are: [Gamma, H_ind, H_omega_lo, H_omega_hi ]
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            M = ModelStructure;
            % Local variables
            p = M.p;
            n_eps = M.n_eps;
            n_zeta = M.n_zeta;
            MP_eps = M.MP_eps;
            MP_zeta = M.MP_zeta;
            n_a_fine = M.n_a_fine;
            a_grid_fine=M.a_grid_fine;
            theta_a = M.theta_a;
            theta_a_f = M.theta_a_f;
            G_ap_fine = M.G_ap_fine;

            a_min = M.p.a_min;
            Hist_max_iter = M.p.Hist_max_iter;
            Hist_tol = M.p.Hist_tol;
            Hist_eta = M.p.Hist_eta;

            

            disp('----------------------------------------------------------------');
            disp('Beginning Histogram Method with Loops');
            disp('----------------------------------------------------------------');

            % Change max iter
            if isempty(N_H)==1
                N_H = Hist_max_iter;
            end

            % Initial distribution
            if isempty(Gamma_0)==1
                Gamma_0 = M.Gamma;
            end

            % Discretize distribution
            disp('Discretizing Choices and Computing Transition Probabilities');
            H_ind = randi([0,1], [n_a_fine,n_eps,n_zeta]);
            H_omega_lo = randi([0,1], [n_a_fine,n_eps,n_zeta,n_eps,n_zeta]);
            H_omega_hi = randi([0,1], [n_a_fine,n_eps,n_zeta,n_eps,n_zeta]);
            a_max = max(a_grid_fine);

            for iter_zeta=1:n_zeta % Current zeta
                Pr_zetap = transpose(MP_zeta.Pi(iter_zeta,:)); % Transitions of future zeta conditional on current zeta
            for iter_eps = 1:n_eps % Current epsilon
                Pr_epsp = transpose(MP_eps.Pi(iter_eps,:)); % Transitions of future epsilon conditional on current epsilon
            for iter_a = 1:n_a_fine % Current a

                % Get index and weight of lower bound on approximation interval
                if (G_ap_fine(iter_a,iter_eps,iter_zeta)>= a_max)
                    H_ind(iter_a, iter_eps, iter_zeta) = n_a_fine - 1 ;
                    omega_lo = 0;
                elseif (G_ap_fine(iter_a,iter_eps,iter_zeta)<= a_min)
                    H_ind(iter_a, iter_eps,iter_zeta) = 1;
                    omega_lo = 1;
                else
                    H_ind(iter_a, iter_eps,iter_zeta) = VFI_Toolbox.Grid_Inv(G_ap_fine(iter_a, iter_eps, iter_zeta), n_a_fine, theta_a_f, a_min, a_max, "Poly" );
                    omega_lo = 1-min(1, max(0,(G_ap_fine(iter_a, iter_eps, iter_zeta)- a_grid_fine(H_ind(iter_a, iter_eps, iter_zeta)))/(a_grid_fine(H_ind(iter_a, iter_eps, iter_zeta)+1)-a_grid_fine(H_ind(iter_a, iter_eps, iter_zeta)))));
                end

                % Store weights for lower and upper bounds on approximation
                % interval, including transition to future states 

                for iter_zetap=1:n_zeta % Future zeta
                    for iter_epsp=1:n_eps % Future eps
                        H_omega_lo(iter_a,iter_eps,iter_zeta,iter_epsp, iter_zetap) = (  omega_lo)*Pr_epsp(iter_epsp)*Pr_zetap(iter_zetap);
                        H_omega_hi(iter_a,iter_eps,iter_zeta,iter_epsp, iter_zetap) = (1-omega_lo)*Pr_epsp(iter_epsp)*Pr_zetap(iter_zetap);
                    end
                end
            end
            end
            end
             % % Correct corner solutions above
            % H_weight[H_ind.==n_a_fine] .= 0
            % H_ind[H_ind.==n_a_fine]    .= n_a_fine-1
            % # Check bounds for weights
            % H_weight = min.(1,max.(0,H_weight))

            % Loop for updating histogram
            disp('Iterating on the Distribution')
            H_dist = 1;
            for iter_H = 1:N_H
                % Update histogram
                Gamma = zeros(n_a_fine, n_eps, n_zeta);
                for iter_zeta = 1:n_zeta % Current zeta
                for iter_eps  = 1:n_eps % Current eps
                for iter_a    = 1:n_a_fine % Current a
                    iter_ap = H_ind(iter_a,iter_eps, iter_zeta);
                    Gamma_0_i = Gamma_0(iter_a, iter_eps, iter_zeta);
                    for iter_zetap=1:n_zeta % Future zeta
                    for iter_epsp =1:n_eps % Future eps
                        % Update is the product of probabilities by
                        % independence of F(eps) and F(zeta)
                        Gamma(iter_ap    , iter_epsp, iter_zetap) = Gamma(iter_ap    , iter_epsp, iter_zetap) + H_omega_lo(iter_a, iter_eps, iter_zeta, iter_epsp, iter_zetap)*Gamma_0_i;
                        Gamma(iter_ap +1 , iter_epsp, iter_zetap) = Gamma(iter_ap + 1, iter_epsp, iter_zetap) + H_omega_hi(iter_a, iter_eps, iter_zeta, iter_epsp, iter_zetap)*Gamma_0_i;
                    end
                    end
                end
                end
                end
            % Update distance
            H_dist = max(max(max(abs(Gamma - Gamma_0))));
            %disp("-------------------------------------------")
            %fprintf("H_dist is %d \n", H_dist )
            %disp("-------------------------------------------")
            % Update initial distribution
            Gamma_0 = (1-Hist_eta)*Gamma + Hist_eta*Gamma_0;
            % Report progress
            if  mod(iter_H,50)==0
                    fprintf('Histogram Loop: iter = %d, dist = %.4d, E[a] = %.4d \n',iter_H,H_dist, sum(a_grid_fine*squeeze(sum(Gamma,[3,2]) )) )
            end
            % Check convergence
            if (H_dist < Hist_tol)
                fprintf(' Histogram iteration converged in iteration %d. H_dist=%.4d \n',iter_H, H_dist)
                M.Gamma = Gamma;
                M.H_ind = H_ind;
                M.H_omega_lo = H_omega_lo;
                M.H_omega_hi = H_omega_hi;

                ModelStructure_Output = M;
                return
            end
            end % End for iter_H


        end % End of Histogram_Method_Loop(ModelStructure, N_H, Gamma_0)
    
            

% ------------------------------------------------------------------------------------------- 
% -------------------------------------------------------------------------------------------

% Histogram iteration method

function Output1 = Histogram_Iteration(ModelStructure, N_H, Gamma_0) 
    % Output 1 is Gamma_0

        M = ModelStructure;
            % Local variables
        n_eps      = M.n_eps;
        n_zeta     = M.n_zeta;
        n_a_fine   = M.n_a_fine;
        H_ind      = M.H_ind;
        H_omega_lo = M.H_omega_lo;
        H_omega_hi = M.H_omega_hi;

        for iter_H = 1:N_H
            % Update histogram - Loop only though points with enough mass
            Gamma = zeros(n_a_fine, n_eps, n_zeta);
            ind = find(Gamma_0 >=1e-12);
            [row,col,dim] = ind2sub([size(Gamma_0,1) size(Gamma_0,2) size(Gamma_0,3)], ind );
            for iter_ind=1:numel(ind)
                iter_ap = H_ind(row(iter_ind),col(iter_ind),dim(iter_ind) );
                for iter_zetap=1:n_zeta % Future zeta
                for iter_epsp=1:n_eps % Future epsilon
                    Gamma(iter_ap  , iter_epsp, iter_zetap)   = Gamma(iter_ap  , iter_epsp, iter_zetap) + H_omega_lo(row(iter_ind),col(iter_ind),dim(iter_ind), iter_epsp, iter_zetap)*Gamma_0(row(iter_ind),col(iter_ind),dim(iter_ind));
                    Gamma(iter_ap+1, iter_epsp, iter_zetap)   = Gamma(iter_ap+1, iter_epsp, iter_zetap) + H_omega_hi(row(iter_ind),col(iter_ind),dim(iter_ind), iter_epsp, iter_zetap)*Gamma_0(row(iter_ind),col(iter_ind),dim(iter_ind));
                end
                end
            end
            Gamma_0 = Gamma;
            % fprintf('Iteration %d \n', iter_H) 
        end
        %fprintf('End of Histogram Iteration')
        Output1 = Gamma_0;

end % End of function Histogram Iteration


% ------------------------------------------------------------------------------------------- 
% -------------------------------------------------------------------------------------------

% Euler Error Function

% G_ap interpolation
function itp_a = G_ap_epsa(iter_zeta, iter_eps, a, ModelStructure) 
    % a (input) is an array of size n, where the values will be evaluated. 
    M = ModelStructure;
    x_grid = M.a_grid;
    f_grid = M.G_ap(:,iter_eps, iter_zeta);
    f_grid = transpose(f_grid);
    itp_a = spline(x_grid, f_grid, a);

end % End of function G_ap_epsa


% Euler Equation Percentage Error

function Output1 = Euler_Error(iter_eps, iter_zeta, a, ModelStructure)
    
    M = ModelStructure;
    % Local variables
    p = M.p;
    MP_eps = M.MP_eps;
    MP_zeta = M.MP_zeta;
    n_eps = M.n_eps;
    n_zeta = M.n_zeta;
    eps_grid = M.eps_grid;
    zeta_grid = M.zeta_grid;
    beta = M.p.beta;
    r = M.p.r;
    w = M.p.w;
    a_grid = M.a_grid;

    %  Interpolate G_ap at current epsilon
    ap = min(a_grid(end),max(a_grid(1), Functions_ModelSolution.G_ap_epsa(iter_zeta, iter_eps, a, M)));

    % Current Consumption
    c = (1+r*zeta_grid(iter_zeta))*a + w*eps_grid(iter_eps) - ap ;

    % Compute LHS of Euler Equation
    LHS = Functions_ModelSolution.d_utility(c, p);

    % Compute right hand side of Euler equation
    RHS_prob = zeros(n_eps, n_zeta);
        % Marginal utility at epsilon', zeta', a', G_ap(a', epsilon',
        % zeta')
        for iter_epsp  = 1:n_eps
        for iter_zetap = 1:n_zeta
            cp = (1+r*zeta_grid(iter_zetap))*ap + w*eps_grid(iter_epsp) - Functions_ModelSolution.G_ap_epsa(iter_zetap, iter_epsp, a, M);
            up = Functions_ModelSolution.d_utility(cp, p);
            RHS_prob(iter_epsp,iter_zetap) = beta*(MP_eps.Pi(iter_eps, iter_epsp)*MP_zeta.Pi(iter_zeta, iter_zetap))*((1+r*zeta_grid(iter_zetap))*up);
        end
        end
     RHS = sum(RHS_prob, 'all');
     % Return percentage error in Euler equation

    Output1 = (RHS/LHS - 1)*100;

    end % End of function Euler_Error




% ------------------------------------------------------------------------------------------- 
% -------------------------------------------------------------------------------------------

    % Aiyagari Equilibrium

    function ModelStructure_Output = Aiyagari_Equilibrium(ModelStructure_in)
    
        % Input is a model structure
        M_in = ModelStructure_in;
    
        % Read files if required
        if M_in.read_flag == true
            G_ap_0   = reshape( readmatrix('File_Folder/Policy_Function.csv') , [M_in.n_a, M_in.n_eps, M_in.n_zeta]);
            %Gamma_0 = reshape( readmatrix('File_Folder/Distribution.csv')    , [M_in.n_a_fine, M.n_eps, M.n_zeta] );
        else
            G_ap_0 = zeros(0,0,1);
            % Gamma_0 = zeros(0,0,1);
        end
    
        Gamma_0 = zeros(0,0,1);
    
        % Compute Policy Functions 
            % PFI_Fixed_Point
            M_in = Functions_ModelSolution.PFI_Fixed_Point('T_EGM_G',M_in,G_ap_0);
    
        % Compute Distribution
            % Histogram_Method_Loop 
            N_H = zeros(0,0,1); % This is to obtain an input equivalent to nothing in the Julia code
            M_in  = Functions_ModelSolution.Histogram_Method_Loop(M_in, N_H, Gamma_0);
    
        % Save Results
        writematrix(M_in.G_ap,'./File_Folder/Policy_Function.csv')
        writematrix(M_in.Gamma,'./File_Folder/Distribution.csv')
    
        % ModelStructure output
        ModelStructure_Output = M_in ;

    end % End of function Aiyagari_Equilibrium





         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end % End of Static Methods

end % End of classdef Functions_ModelSolution

