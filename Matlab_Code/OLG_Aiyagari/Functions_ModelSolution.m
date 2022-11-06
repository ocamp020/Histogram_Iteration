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
            
        % Policy Function Iteration: PFI Fixed Point
        function ModelStructure_Output = Solve_Policy_Functions(BellmanOperator,ModelStructure) % 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % BellmanOperator is a string that tells us which function we
            % are calling
            % ModelStructure is the input model structure
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Local variables
            p = ModelStructure.p ;
            n_eps = ModelStructure.n_eps ;
            n_a  = ModelStructure.n_a  ;
            n_a_fine = ModelStructure.n_a_fine ;
            theta_a = ModelStructure.theta_a ;
            a_grid = ModelStructure.a_grid ;
            a_grid_fine = ModelStructure.a_grid_fine ;
            
            r = p.r ;
            %w = p.w; % Not used but I didnt want to remove it just in case.
            Max_Age = p.Max_Age;

            disp(" ")
            disp('-----------------------------')
            fprintf('Backwards Induction - n_eps=%i, Max_Age=%i, n_a=%i - theta_a=%8.4f - r=%8.4f  \n',n_eps, Max_Age, n_a, theta_a, r )
                % Update value function

                if BellmanOperator=="T_BI_G"
                    [G_ap, G_c] = Functions_ModelSolution.T_BI_G(ModelStructure);
                end

                % Check Borrowing Constraint
                if any(G_ap<p.a_min)
                    error("Borrowing Constraint Violated")
                end

                % Interpolate to fine grid
                G_ap_fine = zeros(n_a_fine, n_eps, Max_Age) ;
                G_c_fine  = zeros(n_a_fine, n_eps, Max_Age) ;

                for age   = 1:Max_Age
                for i_eps = 1:n_eps
                    x=a_grid;
                    y = G_ap(:,i_eps,age);
                    sp = fn2fm(spline(x,y),'B-');
                    yy = fnval(sp, a_grid_fine);
                    G_ap_fine(:,i_eps,age) = transpose(yy);
                    clear y sp yy 
                    y = G_c(:,i_eps,age);
                    sp = fn2fm(spline(x,y),'B-');
                    yy = fnval(sp, a_grid_fine);
                    G_c_fine(:,i_eps,age) = transpose(yy);
                    clear x y sp yy 
                end
                end

                % Update Model
                ModelStructure.G_ap = G_ap;
                ModelStructure.G_ap_fine = G_ap_fine;
                ModelStructure.G_c = G_c;
                ModelStructure.G_c_fine = G_c_fine;
                disp("   ")
                disp("------------------------")
                disp("   ")
                ModelStructure_Output = ModelStructure;
                return
        end % End of function Solve_Policy_Functions

% ------------------------------------------------------------------------------------------- 
% -------------------------------------------------------------------------------------------
        
        % Bellman operator - EGM - Backwards Induction


        function [Output1, Output2] = T_BI_G(Model_Structure) 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % When calling this function it should be like:
            % (M.G_ap, M.G_c) = T_BI_G(M) 
            % Where M is the name of the Model Structure.
            % Output1 is G_ap
            % Output2 is G_c
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % p, n_ζ, n_eps, MP_eps, MP_ζ, n_a, G_ap, method = M
            % β, a_min, r, w = p
            
            M = Model_Structure;
            % Local variables
            p = M.p;
            n_eps = M.n_eps;
            MP_eps = M.MP_eps;
            n_a = M.n_a;
            y_mat = M.y_mat;
            a_grid=M.a_grid;
            a_mat_aeps =M.a_mat_aeps ; 
            G_ap = M.G_ap;
            a_mat = M.a_mat;

            beta = p.beta;
            a_min = p.a_min;
            r = p.r;
            w = p.w;
            Max_Age = p.Max_Age;
            Surv_Pr = p.Surv_Pr;
            Age_Pi = p.Age_Pi;

            

            % Allocate policy functions
            G_c = zeros(n_a, n_eps, Max_Age);
            G_ap = zeros(n_a, n_eps, Max_Age);
            % Final Period savings are zero
            G_c(:,:,Max_Age) = (1+r).*a_mat_aeps + y_mat(:,:,Max_Age) - a_min;
            G_ap(:,:,Max_Age) = a_min;

            % Backward Induction Loop
            for age=(Max_Age-1):-1:1
                % Define RHS of Euler equation for each (a', epsilon, age)
                % The matrix d_utility is (a',epsilon') for a given age, the transition matrix is transposed so that it gives (epsilon',epsilon), result is (a',epsilon)
                Euler_RHS = beta*(Surv_Pr(age)*(1+r).* Functions_ModelSolution.d_utility( (1+r).*a_mat_aeps + y_mat(:,:,age+1) - G_ap(:,:,age+1), p )*transpose(MP_eps.Pi) ) ;
                % Check Monotonicity
                if any(Euler_RHS<0)
                    error("RHS must be monotone for EGM to work")
                end
                % Endogenous consumption for all levels of bequests (ap)
                C_endo = Functions_ModelSolution.d_utility_inv(Euler_RHS,p);
                % Define endogenous grid on assets
                A_endo = (C_endo + a_mat_aeps - y_mat(:,:,age))./(1+r);
                % Interpolate functions on exogenous grid
                for i_eps = 1:n_eps
                    % Sort A_endo for interpolation
                    [A_aux,sort_ind] = sort(A_endo(:,i_eps));
                    C_aux = C_endo(:,i_eps);
                    C_aux = C_aux(sort_ind);
                    % Check boundary condition
                    if min(A_aux)>a_min
                        a_vec = transpose(a_grid(a_grid<min(A_aux)));
                        A_aux = [a_vec ; A_aux];
                        C_aux = [((1+r).*a_vec) + y_mat(1,i_eps,age) - a_min ;C_aux];
                    end
                    C_ip = griddedInterpolant(A_aux, C_aux,'spline','nearest');
                    G_c(:,i_eps,age) = C_ip(a_grid);
                    G_ap(:,i_eps,age) = transpose((1+r).*a_grid) + y_mat(:,i_eps,age) - G_c(:,i_eps,age);
                end
            end
            % Adjust for numerical error
            ind = find(abs(G_ap - a_min)<= 1e-10 );
            G_ap(ind) = a_min;
            G_c(ind) = (1+r)*a_mat(ind) + y_mat(ind) - a_min;
            % Check for borrowing constraint
            % if any(G_ap<a_min)
                %error('Borrowing Constraint Violated')
            % end 
            % Return Results
            Output1 = G_ap;
            Output2 = G_c;
            clear M
        end % End of T_BI_G function

 
% ------------------------------------------------------------------------------------------- 
% -------------------------------------------------------------------------------------------
 
% Histogram method loop
% This function updates 

        function ModelStructure_Output = Histogram_Method_Loop(ModelStructure,Gamma_0)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % When calling this function it should be like:
            % ModelStructure_Output = Histogram_Method_Loop(M, Gamma_0) 
            % The outputs of this function are updated values of
            % Gamma, H_ind, H_omega_lo_s, H_omega_hi_s, H_omega_lo_d, H_omega_hi_d 
            % in the ModelStructure_Output structure
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            M = ModelStructure;
            % Local variables
            p = M.p;
            n_eps = M.n_eps;
            MP_eps = M.MP_eps;
            n_a_fine = M.n_a_fine;
            a_grid_fine=M.a_grid_fine;
            theta_a = M.theta_a;
            theta_a_f = M.theta_a_f;
            G_ap_fine = M.G_ap_fine;
            a_max = M.a_max;

            a_min = p.a_min;
            Hist_max_iter = p.Hist_max_iter;
            Hist_tol = p.Hist_tol;
            Hist_eta = p.Hist_eta;
            Max_Age = p.Max_Age;
            Surv_Pr = p.Surv_Pr;
            


            disp('----------------------------------------------------------------');
            disp('Beginning Histogram Method with Loops');
            disp('----------------------------------------------------------------');


            % Initial distribution
            if isempty(Gamma_0)==1
                Gamma_0 = M.Gamma;
            end
                
            % Set median epsilon for newborns
            med_eps = int64(round(n_eps/2));

            % Set initial wealth to (close to) $1k 
            ind_inner = find(a_grid_fine>=1,1);
            n_a_fine_t = transpose(1:1:n_a_fine);
            b_ind = n_a_fine_t(ind_inner);
            clear ind_inner n_a_fine_t

            % Discretize distribution
            disp('Discretizing Choices and Computing Transition Probabilities');
            H_ind = randi([0,1], [n_a_fine,n_eps,Max_Age]);
            H_omega_lo_s = randi([0,1], [n_a_fine,n_eps,Max_Age, n_eps]);
            H_omega_hi_s = randi([0,1], [n_a_fine,n_eps,Max_Age, n_eps]);
            H_omega_lo_d = randi([0,1], [n_a_fine,n_eps,Max_Age]);
            H_omega_hi_d = randi([0,1], [n_a_fine,n_eps,Max_Age]);
            
            for age=1:Max_Age % Current Age
                Pr_age = Surv_Pr(age) ; % Transitions of future age conditional on current age

                for i_eps = 1:n_eps % Current epsilon
                    Pr_epsp = transpose(MP_eps.Pi(i_eps,:)); % Transitions of future epsilon conditional on current epsilon

                    for i_a=1:n_a_fine % Current a

                        % Get index and weight of lower bound on approximation interval

                        % Get index and weight of lower bound on approximation interval
                        if (G_ap_fine(i_a,i_eps,age)>= a_max)
                            H_ind(i_a, i_eps, age) = n_a_fine - 1 ;
                            omega = 0;
                        elseif (round(G_ap_fine(i_a,i_eps,age),3)<= a_min)
                            H_ind(i_a, i_eps,age) = 1;
                            omega = 1;
                        else
                            H_ind(i_a, i_eps, age) = VFI_Toolbox.Grid_Inv(G_ap_fine(i_a, i_eps, age), n_a_fine, theta_a_f, a_min, a_max, "Poly" );
                            omega = min(1, max(0,(G_ap_fine(i_a, i_eps, age)- a_grid_fine(H_ind(i_a, i_eps, age)))/(a_grid_fine(H_ind(i_a, i_eps, age)+1)-a_grid_fine(H_ind(i_a, i_eps, age)))));
                        end

                        % Transition matrices for survival
                        for i_epsp = 1:n_eps
                            H_omega_lo_s(i_a,i_eps,age,i_epsp) = Pr_age*Pr_epsp(i_epsp)*omega;
                            H_omega_hi_s(i_a,i_eps,age,i_epsp) = Pr_age*Pr_epsp(i_epsp)*(1-omega);
                        end

                        % Transition matrices for death
                        H_omega_lo_d(i_a,i_eps,age) = (1-Pr_age)*omega;
                        H_omega_hi_d(i_a,i_eps,age) = (1-Pr_age)*(1-omega);

                        % Check for conservation of mass
                        H_res = sum(sum( H_omega_lo_s(i_a,i_eps,age,:)+H_omega_hi_s(i_a,i_eps,age,:) ) + H_omega_lo_d(i_a,i_eps,age)) ;
                        if abs((1-H_res) - H_omega_hi_d(i_a,i_eps,age))<1e-5
                            H_omega_hi_d(i_a,i_eps,age) = 1 - H_res;
                        else
                            fprintf("Error in transition probabilities: i_a=%i, i_eps=%i, age=%i, residual=%8.4f \n", i_a, i_eps, age, abs( (1-H_res) - H_omega_hi_d(i_a,i_eps,age)))
                        end
                    end % End loop on assets
                end % End loop on epsilon
            end % End loop on age 
            
            % Check transitions
            H_aux = sum(H_omega_lo_s + H_omega_hi_s, 4)+ H_omega_lo_d + H_omega_hi_d ;
            fprintf(" Check transition functions: maximum = %4.4f, minimum = %4.4f \n", max(max(max(H_aux)-1)), min(min(min(H_aux)-1)))
            H_omega_hi_d = 1 - (sum(H_omega_lo_s+H_omega_hi_s,4)+H_omega_lo_d);
            fprintf(" Check initial distribution: sum(Gamma_0) = %4.4f  \n ", sum(sum(sum(Gamma_0))))

            % Loop for updating histogram
            disp("         ")
            disp(" Iterating on the Distribution ")
            H_dist=1;
            for i_H=1:Hist_max_iter
                % Update histogram
                Gamma = zeros(n_a_fine, n_eps, Max_Age);

                % Update until second to last period of life
                for age=1:Max_Age-1 % Current age
                    for i_eps=1:n_eps % Current epsilon
                        for i_a=1:n_a_fine % Current a

                            i_ap=H_ind(i_a,i_eps, age);

                            % If agents survive
                            for i_epsp=1:n_eps % Future eps
                                Gamma(i_ap  ,i_epsp,age+1) = Gamma(i_ap  ,i_epsp,age+1) + H_omega_lo_s(i_a,i_eps,age,i_epsp)*Gamma_0(i_a,i_eps,age);
                                Gamma(i_ap+1,i_epsp,age+1) = Gamma(i_ap+1,i_epsp,age+1) + H_omega_hi_s(i_a,i_eps,age,i_epsp)*Gamma_0(i_a,i_eps,age);
                            end

                            % If agents die: age=1 and eps=median(eps) and no wealth 
                            Gamma(b_ind, med_eps, 1) = Gamma(b_ind, med_eps,1) + (1-Surv_Pr(age))*Gamma_0(i_a,i_eps,age);
                        end
                    end
                end

                % Final period of life (all agents die)
                for i_eps = 1:n_eps % Current eps
                    for i_a=1:n_a_fine % Current a

                        Gamma(b_ind, med_eps, 1) = Gamma(b_ind, med_eps,1) + Gamma_0(i_a, i_eps,age);
                    end
                end
                % Update distance
                H_dist = max(max(max(abs(Gamma-Gamma_0))));
                
                %Update  initial distribution
                Gamma_0 = (1-Hist_eta)*Gamma + Hist_eta*Gamma_0;

                % Report progress
                
                if  mod(i_H,10)==0
                        fprintf('Histogram Loop: iter = %i, dist = %.2d,  check Gamma=%.2f , check E(a) = %.2f \n',i_H,H_dist , sum(sum(sum(Gamma))) , sum(a_grid_fine*squeeze(sum(Gamma,[3,2]) )) )
                end
                % Check convergence
                if (H_dist < Hist_tol)
                    fprintf(' Histogram iteration converged in iteration %d. H_dist=%.4d \n',i_H, H_dist)
                    M.Gamma = Gamma;
                    M.H_ind = H_ind;
                    M.H_omega_lo_s = H_omega_lo_s;
                    M.H_omega_hi_s = H_omega_hi_s;
                    M.H_omega_lo_d = H_omega_lo_d;
                    M.H_omega_hi_d = H_omega_hi_d;
    
                    ModelStructure_Output = M;
                    clear M
                    return
                end
            end
        end % End of Function Histogram_Method_Loop
           

% ------------------------------------------------------------------------------------------- 
% -------------------------------------------------------------------------------------------

% Histogram iteration method

function Output1 = Histogram_Iteration(ModelStructure, N_H, Gamma_0) 
    % Output 1 is Gamma_0

        
        % Local variables
        M = ModelStructure;
        p=M.p;
        n_eps=M.n_eps;
        n_a_fine = M.n_a_fine; 
        H_ind = M.H_ind; 
        H_omega_lo_s=M.H_omega_lo_s; 
        H_omega_hi_s = M.H_omega_hi_s ;
        a_grid_fine = M.a_grid_fine;
        n_a_fine = M.n_a_fine;

        Max_Age=p.Max_Age;
        Surv_Pr=p.Surv_Pr;
        Gamma_0_inner = Gamma_0;

        % Set median epsilon for newborns
        med_eps = int64(round(n_eps/2));

        % Set initial wealth to (close to) $1k 
        ind_inner = find(a_grid_fine>=1,1);
        n_a_fine_t = transpose(1:1:n_a_fine);
        b_ind = n_a_fine_t(ind_inner);
        clear ind_inner n_a_fine_t

        for iter_H = 1:N_H
            % Update histogram - Loop only though points with enough mass
            Gamma_HI = zeros(n_a_fine, n_eps, Max_Age);
            ind = find(Gamma_0_inner >=1e-12);
            [row,col,ndim] = ind2sub([size(Gamma_0_inner,1) size(Gamma_0_inner,2) size(Gamma_0_inner,3)], ind );
            for iter_ind=1:numel(ind)
                if ndim(iter_ind)<Max_Age % Before Max Age
                % Savings at current state
                i_ap = H_ind(row(iter_ind),col(iter_ind),ndim(iter_ind) );
                
                % If agents survive
                for i_epsp=1:n_eps % Future epsilon
                    Gamma_HI(i_ap  ,i_epsp,ndim(iter_ind)+1) = Gamma_HI(i_ap  ,i_epsp,ndim(iter_ind)+1) + H_omega_lo_s(row(iter_ind), col(iter_ind), ndim(iter_ind),i_epsp)*Gamma_0_inner(row(iter_ind), col(iter_ind), ndim(iter_ind));
                    Gamma_HI(i_ap+1,i_epsp,ndim(iter_ind)+1) = Gamma_HI(i_ap+1,i_epsp,ndim(iter_ind)+1) + H_omega_hi_s(row(iter_ind), col(iter_ind), ndim(iter_ind),i_epsp)*Gamma_0_inner(row(iter_ind), col(iter_ind), ndim(iter_ind));
                end

                %  If agents die: age=1 and epsilon=median(epsilon) and no wealth  
                Gamma_HI(b_ind,med_eps,1) = Gamma_HI(b_ind,med_eps,1) + (1-Surv_Pr(ndim(iter_ind)))*Gamma_0_inner(row(iter_ind), col(iter_ind), ndim(iter_ind));

                else % Max Age, everyone dies 
                    Gamma_HI(b_ind,med_eps,1) = Gamma_HI(b_ind,med_eps,1) + Gamma_0_inner(row(iter_ind), col(iter_ind), ndim(iter_ind));
                end
            end
            % Update distribution
            Gamma_0_inner = Gamma_HI;
            % fprintf('Iteration %d \n', iter_H) 
        end
        %fprintf('End of Histogram Iteration')
        Output1 = Gamma_0_inner;

end % End of function Histogram Iteration


% ------------------------------------------------------------------------------------------- 
% -------------------------------------------------------------------------------------------

% Euler Error Function

% G_ap interpolation
function itp_a = G_ap_epsa(age, iter_eps, a, ModelStructure) 
    % a (input) is an array of size n, where the values will be evaluated. 
    M = ModelStructure;
    x_grid = M.a_grid;
    f_grid = M.G_ap(:,iter_eps, age);
    f_grid = transpose(f_grid);
    itp_a = spline(x_grid, f_grid, a);

end % End of function G_ap_epsa


% Euler Equation Percentage Error

function Output1 = Euler_Error(iter_eps, age, a, ModelStructure)
    
    M = ModelStructure;
    % Local variables
    p = M.p;
    MP_eps = M.MP_eps;
    n_eps = M.n_eps;
    eps_grid = M.eps_grid;
    y_mat = M.y_mat;
    a_grid = M.a_grid;


    beta = p.beta;
    r = p.r;
    w = p.w;
    Max_Age=p.Max_Age;
    Surv_Pr = p.Surv_Pr;


    %  Interpolate G_ap at current epsilon
    ap = min(a_grid(end),max(a_grid(1), Functions_ModelSolution.G_ap_epsa(age, iter_eps, a, M)));

    % Current Consumption
    c = (1+r)*a + y_mat(1,iter_eps,age) - ap ;

    % Compute LHS of Euler Equation
    LHS = Functions_ModelSolution.d_utility(c, p);

    % Compute right hand side of Euler equation
    if age<Max_Age
        RHS_prob = zeros(n_eps,1);
        % Marginal utility at age+1, epsilon', a', G_ap(a', epsilon',age+1)
        for iter_epsp  = 1:n_eps
            cp = (1+r)*ap + y_mat(1,iter_epsp,age+1) - Functions_ModelSolution.G_ap_epsa(age+1, iter_epsp, ap, M);
            up = Functions_ModelSolution.d_utility(cp, p);
            RHS_prob(iter_epsp) = beta*(Surv_Pr(age)*MP_eps.Pi(i_eps,i_epsp)*(1+r)*up);
        end
        RHS = sum(RHS_prob, 'all');
    else
        RHS=0.0;
    end
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
            Gamma_0 = reshape( readmatrix('File_Folder/Distribution.csv')    , [M_in.n_a_fine, M.n_eps, M.P.Max_Age] );
        else
            Gamma_0 = zeros(0,0,1);
        end
    
        % Compute Policy Functions 
            % PFI_Fixed_Point
            M_in = Functions_ModelSolution.Solve_Policy_Functions('T_BI_G',M_in);
        
        for j=1:M_in.p.Max_Age
            M_in.G_ap(:,:,j) = readmatrix( sprintf(['/Users/cyberdim/Dropbox/' ...
                'WESTERN_ECONOMICS/RA_Baxter/Matlab_Hist_Iter/' ...
                'OLG/File_Folder/G_ap_Julia%i.csv'],j));
            M_in.G_ap_fine(:,:,j) = readmatrix( sprintf(['/Users/cyberdim/Dropbox/' ...
                'WESTERN_ECONOMICS/RA_Baxter/Matlab_Hist_Iter/' ...
                'OLG/File_Folder/G_ap_fine_Julia%i.csv'],j));
            M_in.G_c(:,:,j) = readmatrix( sprintf(['/Users/cyberdim/Dropbox/' ...
                'WESTERN_ECONOMICS/RA_Baxter/Matlab_Hist_Iter/' ...
                'OLG/File_Folder/G_c_Julia%i.csv'],j));
            M_in.G_c_fine(:,:,j) = readmatrix( sprintf(['/Users/cyberdim/Dropbox/' ...
                'WESTERN_ECONOMICS/RA_Baxter/Matlab_Hist_Iter/' ...
                'OLG/File_Folder/G_ap_fine_Julia%i.csv'],j));

        end


        % Compute Distribution
            % Histogram_Method_Loop 
            M_in  = Functions_ModelSolution.Histogram_Method_Loop(M_in, Gamma_0);

            for j=1:M_in.p.Max_Age
            M_in.Gamma(:,:,j) = readmatrix( sprintf(['/Users/cyberdim/Dropbox/' ...
                'WESTERN_ECONOMICS/RA_Baxter/Matlab_Hist_Iter/' ...
                'OLG/File_Folder/Gamma_Julia%i.csv'],j));

            end
        % Save Results
        %writematrix(M_in.G_ap,'./File_Folder/Policy_Function.csv')
        %writematrix(M_in.Gamma,'./File_Folder/Distribution.csv')
    
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

