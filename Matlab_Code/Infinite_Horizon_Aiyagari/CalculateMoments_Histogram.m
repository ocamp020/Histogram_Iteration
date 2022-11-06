%###################################################################
%###################################################################
%###################################################################
%## Compute moments using the histogram iteration method 


%###################################################################
%###################################################################
%## 10 Year transition rates for decile 


fprintf('========================================= \n')
fprintf('Computing transition matrix for deciles \n')

% Set time horizon
n_H = 9; % This goes n_H+1 periods ahead

% Find deciles
deciles_a = zeros(11,1);
deciles_a(1) = 0;
for i=2:10
    index = find(100*CDF_a>=(10*(i-1)));
    deciles_a(i) = index(1);
end
deciles_a(11) = M_Aiyagari.n_a_fine;

% For each decile compute transitions   
Tr_deciles_a = zeros(10,10);
for i=1:10
    fprintf('Decile: %d \n', i)
    % Conditional distribution on current decile 
    Gamma_0 = zeros(M_Aiyagari.n_a_fine, M_Aiyagari.n_eps, M_Aiyagari.n_zeta);
    Gamma_0(deciles_a(i)+1:deciles_a(i+1), :, :) = M_Aiyagari.Gamma(deciles_a(i)+1:deciles_a(i+1), :, :);
    Gamma_0 = Gamma_0./sum(sum(sum(Gamma_0)));
    % Iterate distribution
    Gamma_N = Functions_ModelSolution.Histogram_Iteration(M_Aiyagari, n_H, Gamma_0);
    % Obtain transitions (sum over new distribution within deciles)
    inner_column = zeros(10,1);
    for j=1:10
        inner_column(j) = sum(sum(sum(Gamma_N(deciles_a(j)+1: deciles_a(j+1),:,:))));
    end
    Tr_deciles_a(i,:) = inner_column;
    clear inner_column;
end

% Print transition matrix
fprintf(' Wealth Deciles Transition Matrix (N=%d) \n ', n_H)
for i=1:10
    v1 = round(M_Aiyagari.a_grid_fine(deciles_a(i+1)),2);
    v2 = round(100*Tr_deciles_a(i,:),2);
    fprintf('Decile %i, a<=$ %4.2f k // Transition: [%4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f] \n',i,v1,v2)
end

clear v1 v2;

% Save transition matrix to file

fprintf('===============================================   \n')


%###################################################################
%###################################################################
%## 5 Year autocorrelation of consumption

fprintf('=============================================== \n')
fprintf('Autocorrelation of Consumption: 1st quintile    \n')


% Set time horizon
n_H = 2; %This goes n_H+1 periods ahead

% Get average consumption and standard deviation
av_c = sum(sum(sum((M_Aiyagari.G_c_fine.*M_Aiyagari.Gamma))));
sd_c = sqrt(sum(sum(sum( ((M_Aiyagari.G_c_fine-av_c).^2).*M_Aiyagari.Gamma  ))));
% First quintile distribution
Gamma_q = zeros(M_Aiyagari.n_a_fine, M_Aiyagari.n_eps, M_Aiyagari.n_zeta);
Gamma_q(1:deciles_a(3),:,:) = M_Aiyagari.Gamma(1:deciles_a(3),:,:);
Gamma_q = Gamma_q/sum(sum(sum(Gamma_q)));
% First quintile consumption
av_c_q = sum(sum(sum( M_Aiyagari.G_c_fine(1:deciles_a(3),:,:).*Gamma_q(1:deciles_a(3),:,:) )));
sd_c_q = sqrt(sum(sum(sum(  ((M_Aiyagari.G_c_fine(1:deciles_a(3),:,:) - av_c_q).^2).*Gamma_q(1:deciles_a(3),:,:)))));
% Future average and standard deviation conditional on first quintile
    % Iterate distribution
    Gamma_qN = Functions_ModelSolution.Histogram_Iteration(M_Aiyagari, n_H, Gamma_q);
    % Obtain consumption
    av_c_q_N = sum(sum(sum( M_Aiyagari.G_c_fine.*Gamma_qN)));
    sd_c_q_N = sqrt(sum(sum(sum(((M_Aiyagari.G_c_fine-av_c_q_N).^2).*Gamma_qN ))));




% -------------------------------------- commented section in Julia Code -------------------------------------- 
% 
% For each decile compute future consumption
% av_c_deciles = zeros(10)
% sd_c_deciles = zeros(10)
% for i=1:10
    %av_c_deciles(1)= sum(sum(sum(M_Aiyagari.G_c_fine(deciles_a(i)+1:deciles_a(i+1),:,:).*M_Aiyagari.Gamma(deciles_a(i)+1:deciles_a(i+1),:,:))))/sum(sum(sum(M_Aiyagari.Gamma(deciles_a(i)+1:deciles_a(i+1),:,:))));

% ----------------------------------- end of commented section in Julia Code -----------------------------------


% Compute integrand of correlation for first and tenth deciles
cov_c = zeros(M_Aiyagari.n_a_fine, M_Aiyagari.n_eps, M_Aiyagari.n_zeta);
for iter_zeta = 1:M_Aiyagari.n_zeta      % Current zeta
    for iter_eps = 1:M_Aiyagari.n_eps    % Current eps
        % fprintf("   iter_zeta = %i  - iter_eps = %i   \n", iter_zeta, iter_eps);
        % First Quintile
        for iter_a = 1:deciles_a(3) % Current a
            % For each state in state space (a,eps,zeta) get conditional
            % distribution
            Gamma_0 = zeros(M_Aiyagari.n_a_fine, M_Aiyagari.n_eps, M_Aiyagari.n_zeta);
            Gamma_0(iter_a, iter_eps, iter_zeta) = 1 ;
            % Iterate distribution
            Gamma_N = Functions_ModelSolution.Histogram_Iteration(M_Aiyagari,n_H,Gamma_0);
            % Get portion of integrand (C_0 - av_c)*(C_n - av_c)*Gamma_n
            cov_c(iter_a, iter_eps, iter_zeta) = sum(sum(sum( (M_Aiyagari.G_c_fine(iter_a, iter_eps, iter_zeta) - av_c_q ).*(M_Aiyagari.G_c_fine - av_c_q_N).*Gamma_N )));
        end
    end
end

% Integrate with respect to initial distribution of first quintile
cov_c_q = sum(sum(sum(cov_c.*Gamma_q)));
cor_c_q = cov_c_q/sqrt((sd_c_q^2).*(sd_c_q_N^2));

% Print result
disp('                     ')
fprintf('Average Consumption: all - $%8.2fk   \n', round(av_c,2) )
disp('                     ')
fprintf('Average Consumption: q1 - $%8.2fk   \n', round(av_c_q(1),2) )
disp('                     ')
fprintf('Standard Deviation Consumption: all - $%8.2fk   \n', round(sd_c,2) )
disp('                     ')
fprintf('Standard Deviation Consumption: q1 - $%8.2fk   \n', round(sd_c_q(1),2) )
disp('                     ')
fprintf('%i years auto-corr: q1 - $%8.3f  \n', n_H, round(cor_c_q,3))

   
























