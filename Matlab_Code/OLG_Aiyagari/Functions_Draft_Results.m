% Functions
classdef Functions_Draft_Results
    methods(Static)
    
     
% ------------------------------------------------------------------------------------------- 
% -------------------------------------------------------------------------------------------
   % Function H_Moments_Wealth_Profile
   function [pct_a_age_nb, pct_a_age_45] = H_Moments_Wealth_Profile(Model, pct_list, age_0)
       % Inputs are:  Model, pct_list, age_0
       
       % Local variables
       M=Model;
       %%% Define marginal distributions 

 
    Gamma_a = squeeze(sum(M_Aiyagari.Gamma, [3,2])); % Assets 
    Gamma_eps = squeeze(sum(M_Aiyagari.Gamma, [1,3])); % Labor Efficiency 
    Gamma_age = squeeze(sum(M_Aiyagari.Gamma, [1,2])); % Age

%%% Define index for median shocks
    med_eps = int64(round(M_Aiyagari.n_eps/2));

    % Unconditioned Wealth Profile
    pct_a_age_nb = zeros(M.p.Max_Age, length(pct_list)+1);
    for h=1:M.p.Max_Age
        Gamma_aux = sum(M.Gamma(:,:,h), [2,3]);
        Gamma_aux = Gamma_aux(:,1);
        CDF_h = cumsum(Gamma_aux/sum(Gamma_aux));
        for p = 1:numel(pct_list)
            inner = find((100*CDF_h)>=(pct_list(p)) );
             pct_a_age_nb(h,1:end-1) = 
        end
           

            
     end
     clear inner
     pct_a_age_nb(h,1:end) = sum(M.a_mat_fine(:,:,h).*M.Gamma(:,:,h))/Gamma_age(h);

     % 45-year olds with more than median income
     pct_a_age_45 = zeros(M.p.Max_Age-(age_0-1), numel(pct_list)+1);
     % Initial value
     Gamma_h = zeros(size(M.Gamma));
     Gamma_h(:, med_eps, age_0) = M.Gamma(:,med_eps:end,age_0)/sum(M.Gamma(:,med_eps:end,age_0));
     CDF_1 = sum(Gamma_h, [2,3]);
     CDF_1 = CDF_1(:,1);
     CDF_1 = cumsum(CDF_1);
     for p=1:numel(pct_list)
         inner = ;
         pct_a_age_45(1,1:end-1) = 
     end
     clear inner
     pct_a_age_45(1,end) = sum(M.a_mat_fine(:,:,age_0))*Gamma_h(:,:,age_0);
     % Simulate cohort forward
     for h=2:M.p.Max_Age-(age_0-1)
         Gamma_h = Functions_ModelSolution.Histogram_Iteration(M,1,Gamma_h);
         Gamma_h(:,:,(1:end ~= (h+(age_0-1)))) = 0;
         Gamma_h = Gamma_h/sum(Gamma_h);
         CDF_h = sum(Gamma_h, [2,3]);
         CDF_h = CDF_h(:,1);
         CDF_h = cumsum(CDF_h);
         for p=1:numel(pct_list)
             inner = ;
             pct_a_age_45(h,1:end-1)=
         end
         pct_a_age_45(h,end)=sum(M.a_mat_fine(:,:,h)*Gamma_h(:,:,h));
     end
    return

   end % End of function H_Moments_Wealth_Profile
    

     
% ------------------------------------------------------------------------------------------- 
% -------------------------------------------------------------------------------------------
   



% Function H_Moments_Wealth_Corr
function cor_a = H_Moments_Wealth_Corr(M, age_0, age_T)
    % Inputs are M (Model), age_0, age_t

    n_H = age_T - age_0 ;% Simulate for n_H years 
    % Get a_grid vector
    a_grid_vec = M.a_grid_fine;
    % Turn off death to ensure balance panel 
    M_C = M;
    M_C.p.Surv_Pr =  vertcat(ones(M.p.Max_Age -1 , 0));
    % Initial distribution
    Gamma_h = zeros(size(M.Gamma));
    Gamma_h(:,:,age_0) = M.Gamma(:,:,age_0)/sum(M.Gamma(:,:,age_0));
    % Average and standard deviation at initial distribution 
    av_a_0 = sum(a_grid_vec * sum(Gamma_h,[2,3]));
    sd_a_0 = sqrt( sum ((a_grid_vec-av_a_0).^2)*sum(Gamma_h,[2,3]) );
    av_a_N = sum(a_grid_vec * sum(Gamma_h(:,:,n_H+age_0,[2,3])))/sum(M.Gamma(:,:,n_H+age_0));
    sd_a_N = sqrt(sum(((a_grid_vec-av_a_N).^2)*sum(M.Gamma(:,:,n_H+age_0),[2,3])/sum(M.Gamma(:,:,n_H+age_0))));

    % Follow each initial state and fill in integrand
    cov_a = zeros(M.n_a_fine, M.n_eps);
    for i_eps=1.M.n_eps
        for i_a=1:M.n_a_fine
            % For each state in state space(a,eps,h) get conditional
            % distribution.
            Gamma_0 = zeros(M.n_a_fine, M.n_eps, M.p.Max_Age);
            Gamma_0(i_a,i_eps, age_0) = 1;
            % Iterate distribution
            Gamma_N = Functions_ModelSolution.Histogram_Iteration(M_C,n_H,Gamma_0);
            % Fill in integrand (a_0 - av_a)*(a_n - av_a)*Gamma_N
            cov_a(i_a,i_eps) = sum((a_grid_vec(i_a)-av_a_0)*(a_grid_vec-av_a_N)*sum(Gamma_N,[2,3])/sum(Gamma_N));
        end
    end
    % Integrate covariance
    cor_a = sum(cov_a*Gamma_h(:,:,age_0))/sqrt((sd_a_0^2)*(sd_a_N^2));
    % Return moment
    return
end












 
% ------------------------------------------------------------------------------------------- 
% -------------------------------------------------------------------------------------------
   

    end % End of (Static) methods 
end % End of classdef