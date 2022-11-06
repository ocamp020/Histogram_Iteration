%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute moments using the histogram iteration method 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Cohort evolution of 90-10 percentiles of wealth 

disp("===============================================")
disp("90-10 Percentiles for Cohorts")

%%% Newborns
pct_C_age_1 = zeros(2, M_Aiyagari.p.Max_Age);
for h=1:M_Aiyagari.p.Max_Age
    Gamma_aux = sum(M_Aiyagari.Gamma(:,:,h),[2,3]) ;
    CDF_h = cumsum(Gamma_aux/sum(Gamma_aux)) ;
    ind_10 = find((100*CDF_h)>=(10));
    ind_90 = find((100*CDF_h)>=(90));
    pct_C_age_1(1,h) = M_Aiyagari.a_grid_fine(ind_10(1));
    pct_C_age_1(2,h) = M_Aiyagari.a_grid_fine(ind_90(1));
    clear ind_10 ind_90 ;
end

%%% 45-year olds with more than median income
age_0 = 26 ;
pct_C_age_2 = zeros(2, M_Aiyagari.p.Max_Age-(age_0-1)) ;  
% Initial Value
Gamma_h = zeros(size(M_Aiyagari.Gamma));
Gamma_h(:,med_eps:end,age_0) = M_Aiyagari.Gamma(:,med_eps:end,age_0)/(sum(sum(M_Aiyagari.Gamma(:,med_eps:end,age_0))))    ;
CDF_1 = cumsum(sum(Gamma_h,[2,3])) ;
ind_10 = find((100*CDF_1)>=10);
ind_90 = find((100*CDF_1)>=(90));
pct_C_age_2(1,1) = M_Aiyagari.a_grid_fine(ind_10(1));
pct_C_age_2(2,1) = M_Aiyagari.a_grid_fine(ind_90(1));
% Simulate cohort forward
for h=2:(M_Aiyagari.p.Max_Age -(age_0 - 1))
   Gamma_h_inner = Functions_ModelSolution.Histogram_Iteration(M_Aiyagari,1,Gamma_h);
   for h_inner = 1:M_Aiyagari.p.Max_Age
       if h_inner ~=(h+(age_0-1))
           Gamma_h_inner(:,:,h_inner)=0;
       else
           fprintf('success is %i \n', h_inner)
       end
   end

   Gamma_h_inner = Gamma_h_inner / sum(sum(sum(Gamma_h_inner)));

CDF_h = cumsum(sum(Gamma_h_inner,[2,3]));
clear ind_10 ind_90 ;
ind_10 = find((100*CDF_h)>=10);
ind_90 = find((100*CDF_h)>=90);
pct_C_age_2(1,h) = M_Aiyagari.a_grid_fine(ind_10(1));
pct_C_age_2(2,h) = M_Aiyagari.a_grid_fine(ind_90(1));
clear ind_10 ind_90;
Gamma_h = Gamma_h_inner;
end

disp("===============================================")


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Auto-correlation of wealth Ages 35 and 65 (conditional on survival)

disp("===============================================")
disp("Autocorr of Wealth: ages 35-65")

age_0 = 16; % Start agents at age 35
n_H = 65-(35-1) ; % Simulate for n_H years 
% Get a_grid vector
a_grid_vec = M_Aiyagari.a_grid_fine;
% Turn off death to ensure balance panel 
M_Aiyagari_C = M_Aiyagari;
M_Aiyagari_C.p.Surv_Pr = [ones(M_Aiyagari.p.Max_Age - 1,1); 0] ;
% Initial distribution
Gamma_h = zeros(size(M_Aiyagari.Gamma));
Gamma_h(:,:,age_0) = M_Aiyagari.Gamma(:,:,age_0)./sum(sum(M_Aiyagari.Gamma(:,:,age_0))) ;
% Average and standard deviation at initial distribution 
av_a_0 = sum(a_grid_vec*sum(Gamma_h,[2,3]));
sd_a_0 = sqrt( sum( ((a_grid_vec-av_a_0).^2)*sum(Gamma_h,[2,3]) ));
av_a_N = sum(a_grid_vec*sum(M_Aiyagari.Gamma(:,:,n_H+age_0),[2,3] ) )/sum(sum(M_Aiyagari.Gamma(:,:,n_H+age_0)));
sd_a_N = sqrt(sum(transpose((a_grid_vec - av_a_N).^2).*sum(M_Aiyagari.Gamma(:,:,n_H+age_0),[2,3])/sum(sum(M_Aiyagari.Gamma(:,:,n_H+age_0)) ) ));
% Follow each initial state and fill in integrand 
cov_a_3565 = zeros(M_Aiyagari.n_a_fine, M_Aiyagari.n_eps);
for i_eps = 1:M_Aiyagari.n_eps
    fprintf('Iterating with i_eps = %i \n', i_eps)
    for i_a = 1:M_Aiyagari.n_a_fine
        % fprintf('Iterating with i_ϵ=%i , i_a=%i \n',i_eps, i_a)
        % For each state in state space (a,eps,zeta) get conditional distribution 
        Gamma_0 = zeros(M_Aiyagari.n_a_fine, M_Aiyagari.n_eps, M_Aiyagari.p.Max_Age);
        Gamma_0(i_a, i_eps,age_0) = 1;
        % Iterate distribution
        Gamma_N = Functions_ModelSolution.Histogram_Iteration(M_Aiyagari_C, n_H, Gamma_0);
        % Fill in integrand (a_0 - av_a)*(a_n - av_a)*Gamma_N
        cov_a_3565(i_a,i_eps) = sum((a_grid_vec(i_a)-av_a_0)*(a_grid_vec - av_a_N)*sum(Gamma_N,[2,3])/sum(sum(sum(Gamma_N))) );
    end
end
% Integrate covariance
cor_a_3565 = sum(sum(cov_a_3565.*Gamma_h(:,:,age_0)))/sqrt((sd_a_0^2)*(sd_a_N^2));
disp("===============================================")











































