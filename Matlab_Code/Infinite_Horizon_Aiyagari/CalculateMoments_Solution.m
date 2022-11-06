%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Define marginal distributions 

 
    Gamma_a = squeeze(sum(M_Aiyagari.Gamma, [3,2])); % Assets 
    Gamma_eps = squeeze(sum(M_Aiyagari.Gamma, [1,3])); % Labor Efficiency 
    Gamma_zeta = squeeze(sum(M_Aiyagari.Gamma, [1,2])); % Returns  

%%% Define index for median shocks
    med_eps = int64(round(M.n_eps/2));
    med_zeta = int64(round(M.n_zeta/2)); 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Print income and return grids

fprintf('============================================================= \n')

%%% Labor income
fprintf('Income Grid and Probability \n')
fprintf('Node - Value - PDF - Gamma \n')
aux = [ transpose(M_Aiyagari.eps_grid)*p.w  transpose(M_Aiyagari.MP_eps.PDF) transpose(Gamma_eps)];
for i_eps=1:M_Aiyagari.n_eps
    %inner = round(aux(i_eps,:),4);
    fprintf('%d : [%8.4f %8.4f %8.4f] \n', i_eps, round(transpose(aux(i_eps,:)),4))
end

v1 = round(dot(aux(:,1), aux(:,2)   ),3 );
v2 = round(dot(aux(:,1), Gamma_eps  ),3 );
v3 = round(sqrt( dot(((aux(:,1)-p.w).^2),Gamma_eps)),3);

fprintf('Expected value: $%8.4f Thousand Dollars \n', v1)
fprintf('Expected value: $%8.4f Thousand Dollars \n', v2)
fprintf('Standard Deviation: $%8.4f Thousand Dollars \n', v3)

clear v1 v2 v3;

%%% Return
fprintf('Return Grid and Probability  \n')
fprintf('Node - Value - PDF - Γ  \n')
aux = [transpose(M_Aiyagari.zeta_grid)*p.r*100  M_Aiyagari.MP_zeta.PDF  Gamma_zeta];
for i_zeta = 1:M_Aiyagari.n_zeta
    %inner = round(aux(i_zeta,:),4);
    fprintf('%d : [%8.4f %8.4f %8.4f] \n', i_zeta, round(transpose(aux(i_zeta,:)),4))
end

fprintf('       \n')

v1 = round(dot(aux(:,1), aux(:,2)   ),4 );
v2 = round(dot(aux(:,1), Gamma_zeta  ),4 );
v3 = round(100*sqrt( dot(((aux(:,1)/100 -p.r).^2),Gamma_zeta)),4);

fprintf('Expected value: %6.3f %% \n', v1)
fprintf('Expected value: %6.3f %% \n ', v2)
fprintf('Standard Deviation: $%6.3f %%  \n ', v3)

clear v1 v2 v3;

%%% Wealth-Weighted Returns

fprintf('Wealth Weighted Returns  \n')
r_mat = p.r*M_Aiyagari.zeta_mat_fine;
aux_pr = M_Aiyagari.a_mat_fine.*M_Aiyagari.Gamma;
aux_pr = aux_pr./sum(sum(sum(aux_pr)));
av_r_w = 100.*sum(sum(sum(r_mat.*aux_pr)));
sd_r_w = 100.*sum(sum(sum(((r_mat - av_r_w/100).^2).*aux_pr)));
v1 = round(av_r_w,4);
v2 = round(sd_r_w,4);
fprintf('Expected Value: %6.4f %%  \n ', v1)
fprintf('Standard Deviation: %6.4f %%  \n', v2)

clear v1 v2;
fprintf('=============================================== \n')



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
v1 = round(av_a, 3);
v2 = round(sd_a,3);


fprintf('Expected Value: $%8.4fk \n', v1)
fprintf('Standard Deviation: $%8.4fk \n', v2)
fprintf(' Top X%%      Level       Share \n')

clear v1 v2;

for i=1:5
    v1 = round(100-Top_shares(i,1),2);
    v2 = round(Top_shares(i,2),3);
    v3 = round(Top_shares(i,3),2);
    fprintf('%4.2f  $ %8.3fk %4.2f%% \n', v1, v2, v3)
end
clear v1 v2 v3;
disp('===============================================')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Longitudinal Moments

% # Get moments from histogram method
run("CalculateMoments_Histogram.m")

% # Get moments from Monte Carlo simulation 
run("CalculateMoments_MonteCarlo.m")


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot Grid and Fine Grid (Zoom at the bottom and the top)

% Plots are WIP



