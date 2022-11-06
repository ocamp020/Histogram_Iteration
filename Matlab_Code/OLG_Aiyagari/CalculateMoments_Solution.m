%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Define marginal distributions 

 
    Gamma_a = squeeze(sum(M_Aiyagari.Gamma, [3,2])); % Assets 
    Gamma_eps = squeeze(sum(M_Aiyagari.Gamma, [1,3])); % Labor Efficiency 
    Gamma_age = squeeze(sum(M_Aiyagari.Gamma, [1,2])); % Age

%%% Define index for median shocks
    med_eps = int64(round(M_Aiyagari.n_eps/2));
    ref_age = 31;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Print income and return grids

fprintf('============================================================= \n')

%%% Labor income
fprintf('Income Grid and Probability \n')
fprintf('Node - Value - PDF - Gamma \n')
aux = [ transpose(M_Aiyagari.eps_grid)*M_Aiyagari.p.w  transpose(M_Aiyagari.MP_eps.PDF) transpose(Gamma_eps)];
for i_eps=1:M_Aiyagari.n_eps
    %inner = round(aux(i_eps,:),4);
    fprintf('%d : [%8.4f %8.4f %8.4f] \n', i_eps, round(transpose(aux(i_eps,:)),4))
end

av_y = sum(sum(sum(M_Aiyagari.y_mat_fine.*M_Aiyagari.Gamma)));
sd_y = sqrt(sum(sum(sum(((M_Aiyagari.y_mat_fine - av_y).^2).*M_Aiyagari.Gamma ))));
age_profile_y = zeros(M_Aiyagari.p.Max_Age,1);
for i=1:M_Aiyagari.p.Max_Age
    age_profile_y(i)= ( sum(sum(sum(M_Aiyagari.y_mat_fine(:,:,i).*M_Aiyagari.Gamma(:,:,i))/Gamma_age(i))) );
end

v1 = round(av_y,3);
v2 = round(sd_y,3);

fprintf('Expected value: $%8.4f Thousand Dollars \n', v1)
fprintf('Standard Deviation: $%8.4f Thousand Dollars \n', v2)

clear v1 v2 ;


%%% Asset Distribution Stats 

fprintf("=============================================== \n")
fprintf(" Asset Distribution Stats  \n")
av_a = dot( M_Aiyagari.a_grid_fine, Gamma_a);
sd_a = sqrt(dot((M_Aiyagari.a_grid_fine(1:end) - av_a).^2, Gamma_a));
CDF_a = cumsum(Gamma_a);
Lorenz_a = cumsum(transpose(M_Aiyagari.a_grid_fine).*Gamma_a)/av_a;
Top_shares = zeros(5,3);
Top_shares(:,1) =  [90;95;99;99.9;99.99];
for i=1:5
    ind_top = find(CDF_a>=(Top_shares(i,1)/100));
    Top_shares(i,2) = M_Aiyagari.a_grid_fine(ind_top(1));
    Top_shares(i,3) = 100*dot(M_Aiyagari.a_grid_fine(ind_top), Gamma_a(ind_top))/av_a;
end

age_profile_a = zeros(M_Aiyagari.p.Max_Age,1);
age_profile_ap = zeros(M_Aiyagari.p.Max_Age,1);

for i=1:M_Aiyagari.p.Max_Age
    age_profile_a(i)  = sum(sum(M_Aiyagari.a_mat_fine(:,:,i).*M_Aiyagari.Gamma(:,:,i))/Gamma_age(i));
    age_profile_ap(i) = sum(sum(M_Aiyagari.G_ap_fine(:,:,i).*M_Aiyagari.Gamma(:,:,i))/Gamma_age(i));
end
B2A_ratio = sum(age_profile_ap.*Gamma_age.*(1 - M_Aiyagari.p.Surv_Pr))/av_a;

v1 = round(av_a, 3);
v2 = round(sd_a,3);


fprintf('Expected Value: $%8.4fk \n', v1)
fprintf('Standard Deviation: $%8.4fk \n', v2)
clear v1 v2;

fprintf(' Top X%%      Level       Share \n')
for i=1:5
    v1 = round(100-Top_shares(i,1),2);
    v2 = round(Top_shares(i,2),3);
    v3 = round(Top_shares(i,3),2);
    fprintf('%4.2f%%  $%8.3fk %4.2f%% \n', v1, v2, v3)
end
clear v1 v2 v3;
fprintf('Bequest to wealth ratio = %4.2f%% \n', round(100*B2A_ratio,2))

disp('===============================================')

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Income Distribution




    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Asset Distribution


% Plot CDF

% Plot Pareto Tail (Above $1 Million)

% Plot Lorenz Curve

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Age Profiles



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Longitudinal Moments




% Get moments from histogram interation method

% Get moments from histogram method
run('CalculateMoments_Histogram.m')

% Get moments from Monte Carlo simulation
%run("CalculateMoments_MonteCarlo.m")


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Solution Properties: Grids and Policy Functions



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot Grid and Fine Grid (Zoom at the bottom and the top)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot Saving Functions (median labor efficiency and interest rate)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Euler Errors plots


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
