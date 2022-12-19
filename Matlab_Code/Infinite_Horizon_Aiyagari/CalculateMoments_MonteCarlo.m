
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute moments using the Monte-Carlo method 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run Simulations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate structure for the panel using Parameters module

%M_P_structure = struct('N_Panel', 500000 , 'T_Panel', 10, 'T_Simul', 1000, 'N_Min', 1000, 'Simul_tol', 1e-2, 'rng_seed', 3489398); % True definition.
%M_P_structure.a_mat    = zeros([M_P_structure.N_Panel, M_P_structure.T_Panel]);
%M_P_structure.c_mat    =  zeros([M_P_structure.N_Panel, M_P_structure.T_Panel]);
%M_P_structure.eps_mat  = zeros([M_P_structure.N_Panel, M_P_structure.T_Panel]);
%M_P_structure.zeta_mat = zeros([M_P_structure.N_Panel, M_P_structure.T_Panel]);
%M_P_structure.t_vec    = zeros([M_P_structure.N_Panel,1]);
    
M_P_structure = Functions_MonteCarlo.PanelStructure(500000 ,10,1000, 1000, 1e-2,3489398); % True definition.

M_P = Functions_MonteCarlo.Simulate_Panel(M_Aiyagari, M_P_structure, 1) ; 
% Value for Seed_Flag==1 to be true, Seed_Flag==0 to be false.


fig_sample = [1000 ; transpose(5000:5000:M_P.N_Panel)];
fig_N      = length(fig_sample)         ;


% Save Results
writematrix(M_P.a_mat,'File_Folder/S_a_mat.csv')
writematrix(M_P.c_mat,'File_Folder/S_c_mat.csv')
writematrix(M_P.eps_mat,'File_Folder/S_e_mat.csv')
writematrix(M_P.zeta_mat,'File_Folder/S_z_mat.csv')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Average Wealth, Percentiles, and Top Wealth Shares 

pct_list = [90;95;99;99.9;99.99];
av_a_S       = zeros(fig_N,1)   ;
pct_S        = zeros(5,fig_N) ;
Top_Shares_S = zeros(5,fig_N) ; 


for i=1:fig_N
    a_sample = M_P.a_mat(1:fig_sample(i),end) ;
    av_a_S(i) = mean(a_sample);
    pct_S(:,i) = prctile(a_sample, pct_list);
    inner_vec = zeros(5,1);
    for p=1:5
        I = find(a_sample>=pct_S(p,i));
        inner_vec(p) = (100*sum(a_sample(I)))/(fig_sample(i)*av_a_S(i));
    end
    Top_Shares_S(:,i) = inner_vec ;
end

clear p i

disp('                                      ')
disp(' Comparing Top Shares and Percentiles ')
disp('                                      ')

disp('   Top X%    Share_Hist  Share_Simul     ')
for i=1:5
    fprintf('  %8.2f%%       %8.2f%%     %8.2f%%   \n',  ...
        round(100-Top_shares(i,1),2), round(Top_shares(i,3),2), ...
        round(Top_Shares_S(i,end),2) )
end


disp('   Top X%    Level_Hist  Level_Simul     ')
for i=1:5
    fprintf('  %8.2f%%       $%9.3fk     $%9.3fk   \n',  ...
        round(100-Top_shares(i,1),2), round(Top_shares(i,2),3), ...
        round(pct_S(i,end),3) )
end



%% Figures: Top percentiles 1% and 0.1%     

f=figure();
plot( fig_sample/1000 , pct_S(3,:),'-o' ,'color', [0 0.4470 0.7410] );
yline(mean(Top_shares(3,2)), 'color', [0.8500 0.3250 0.0980])
ylim([floor(min(pct_S(3,:)/500))*500 ceil(max(pct_S(4,:)/500))*500])
title('99th and 99.9th percentiles','FontSize', 18)
xlabel('Sample Size: Thousands','FontSize',14)
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
set(gca,'XTickLabelMode','auto')
hold on
plot( fig_sample/1000 , pct_S(4,:), '-o','color', [0 0.4470 0.7410] );
yline(mean(Top_shares(4,2)), 'color', [0.8500 0.3250 0.0980])
hold off
f.Units = 'centimeters';
f.PaperUnits = 'centimeters';
f.PaperSize = f.Position(3:4);
f.PaperSize = f.Position(3:4)+0.2; % Add 0.1 cm of margin in each direction
saveas(gcf,'Fig_Folder/Top_pct_Simul.pdf')



%% Figures: Top Shares 1% and 0.1%      

f=figure();
plot( fig_sample/1000 , Top_Shares_S(3,:),'-o' ,'color', [0 0.4470 0.7410]);
yline(mean(Top_Shares_S(3,3)), 'color', [0.8500 0.3250 0.0980])
ylim([0 ceil(max(Top_Shares_S(3,:)/5))*5])
xlim([1, M_P.N_Panel/1000])
hold on
plot( fig_sample/1000 , Top_Shares_S(4,:),'-o' ,'color', [0 0.4470 0.7410]);
yline(mean(Top_Shares_S(4,3)), 'color', [0.8500 0.3250 0.0980])
hold off
title('Top 1% and 0.1% Shares','FontSize', 18)
xlabel('Sample Size: Thousands','FontSize',14)
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
set(gca,'XTickLabelMode','auto')
f.Units = 'centimeters';
f.PaperUnits = 'centimeters';
f.PaperSize = f.Position(3:4);
f.PaperSize = f.Position(3:4)+0.2; % Add 0.1 cm of margin in each direction
saveas(gcf,'Fig_Folder/Top_Share_Simul.pdf')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lorenz Curve

Lorenz_S = zeros(100,4);
sample_vec = [50000, 100000, 250000, 500000] ;  %delete next line when
%delivering the final code bc my Matlab cant run such a big N_Panel
%sample_vec = [2000, 4000, 10000, 20000] ; 
for i=1:4
    a_sample = M_P.a_mat(1:sample_vec(i),end);
    av_a_aux = sum(a_sample);
    p_aux = prctile(a_sample, 1:100) ;
    Lorenz_S_inner = zeros(100,1) ;
    for p=1:100
        I = find(a_sample<=p_aux(p)) ;
        Lorenz_S_inner(p) = 100*sum(a_sample(I))/av_a_aux ;
    end
    Lorenz_S(:,i) = Lorenz_S_inner ;
end


%% Figure with Lorenz Curve for 10k, 50k, 100k    
f=figure();
plot(1:100, 1:100, 'color', [.7 .7 .7])
title('Lorenz curve','FontSize', 18)
ylim([0 100])
xlim([0 100])
hold on
plot(100*CDF_a, 100*Lorenz_a,'MarkerFaceColor', [0 0.4470 0.7410] )
plot(1:100, Lorenz_S(:,1))
plot(1:100, Lorenz_S(:,2))
plot(1:100, Lorenz_S(:,3))
plot(1:100, Lorenz_S(:,4))
hold off
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
set(gca,'XTickLabelMode','auto')
f.Units = 'centimeters';
f.PaperUnits = 'centimeters';
f.PaperSize = f.Position(3:4);
f.PaperSize = f.Position(3:4)+0.2; % Add 0.1 cm of margin in each direction
saveas(gcf,'Fig_Folder/Distribution_Wealth_Lorenz_Simul.pdf')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pareto Tail
    % 10,000 observations
    Pareto_a_10k = M_P.a_mat(1:10000,end) ; % Select first 10.000 observations 
    I = find( Pareto_a_10k >= 5000) ;
    Pareto_a_10k = sort(Pareto_a_10k(I)) ; % Select and sort observations above $1M  
    Pareto_p_10k = ((numel(Pareto_a_10k):-1:1)/numel(Pareto_a_10k)) ;

    % 50,000 observations
    Pareto_a_50k = M_P.a_mat(1:50000,end) ; % Select first 10.000 observations 
    I = find( Pareto_a_50k >= 5000) ;
    Pareto_a_50k = sort(Pareto_a_50k(I)) ; % Select and sort observations above $1M  
    Pareto_p_50k = ((numel(Pareto_a_50k):-1:1)/numel(Pareto_a_50k)) ;
    
    % 100,000 observations
    Pareto_a_100k = M_P.a_mat(1:100000,end) ; % Select first 10.000 observations 
    I = find( Pareto_a_100k >= 5000) ;
    Pareto_a_100k = sort(Pareto_a_100k(I)) ; % Select and sort observations above $1M  
    Pareto_p_100k = ((numel(Pareto_a_100k):-1:1)/numel(Pareto_a_100k)) ;
    
    % 500,000 observations
    Pareto_a_500k = M_P.a_mat(1:500000,end) ; % Select first 10.000 observations 
    I = find( Pareto_a_500k >= 5000) ;
    Pareto_a_500k = sort(Pareto_a_500k(I)) ; % Select and sort observations above $1M  
    Pareto_p_500k = ((numel(Pareto_a_500k):-1:1)/numel(Pareto_a_500k)) ;
    

    % Figure with all Pareto tails (above $5m)
f=figure();
plot( log(grid_1M(1:end-1)/5000), P_coeff*log(grid_1M(1:end-1)/5000), 'color', [0.8500 0.3250 0.0980])
xlim([log(1) log(ceil(M_Aiyagari.a_grid(end)/5000)*1)]);
%xlim([log(1) log(100000/5000)])
xticks(log([1 2 4 8 16]))
xticklabels({'$5m', '$10m', '$20m', '$40m', '$80m'} )  
hold on
scatter(log(grid_1M(1:end-1)/5000)  , log(CCDF_1M(1:end-1)) , 10,'filled', 'MarkerFaceColor', [0 0.4470 0.7410])
scatter(log(Pareto_a_10k(1:end-1)/5000)  , log(Pareto_p_10k(1:end-1)) , 10,'filled')
scatter(log(Pareto_a_50k(1:end-1)/5000)  , log(Pareto_p_50k(1:end-1)) , 10,'filled')
scatter(log(Pareto_a_100k(1:end-1)/5000)  , log(Pareto_p_100k(1:end-1)) , 10,'filled')
scatter(log(Pareto_a_500k(1:end-1)/5000)  , log(Pareto_p_500k(1:end-1)) , 10,'filled')
hold off
title('Distribution Tail','FontSize', 18)
xlabel('Log Assets','FontSize',14)
ylim([ floor(log(CCDF_1M(end-1))/4)*4  0  ])
f.Units = 'centimeters';
f.PaperUnits = 'centimeters';
f.PaperSize = f.Position(3:4);
f.PaperSize = f.Position(3:4)+0.2; % Add 0.1 cm of margin in each direction
saveas(gcf,'Fig_Folder/Distribution_Wealth_Pareto_Simul.pdf')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 10 Year transition rates for decile 
    deciles_a_S = zeros(11,fig_N) ;
    Tr_deciles_a_S = zeros(10,10,fig_N) ;
    for i=1:fig_N
        a_aux_0 = M_P.a_mat(1:fig_sample(i),1)  ;
        a_aux_T = M_P.a_mat(1:fig_sample(i),end) ;
        deciles_a_S(:,i) = prctile(a_aux_T, 0:10:100) ; %  Deciles based on end of sample (hoping for stationariety)
        for p=1:10
            a=deciles_a_S(p);
            b=deciles_a_S(p+1);
            ind_d = find((a <= a_aux_0) & (a_aux_0 <=b)) ;
            a_T = a_aux_T(ind_d); %  Follow them to t=10
            Tr_deciles_a_S_inner = zeros(1,10);
            for pT=1:10
                I2 = (deciles_a_S(pT)<= a_T) & (a_T <= deciles_a_S(pT+1)) ;
                Tr_deciles_a_S_inner(pT) = 100*sum(I2)/numel(ind_d);
            end
            Tr_deciles_a_S(p,:,i) = Tr_deciles_a_S_inner;
            clear Tr_deciles_a_S_inner
        end
    end


% Figure with transitions of bottom decile
f=figure();
plot( fig_sample/1000, reshape(Tr_deciles_a_S(1,:,:),10,101))
yline(100*(Tr_deciles_a(1,:)), 'color', [0.8500 0.3250 0.0980])
ylim([0 ceil(max(max(reshape(Tr_deciles_a_S(1,:,:),10,101))))])
xlim([1 M_P.N_Panel/1000])
title('1st Decile Transitions','FontSize', 18)
xlabel('Sample Size: Thousands','FontSize',14)
f.Units = 'centimeters';
f.PaperUnits = 'centimeters';
f.PaperSize = f.Position(3:4);
f.PaperSize = f.Position(3:4)+0.2; % Add 0.1 cm of margin in each direction
saveas(gcf,'Fig_Folder/Tr_Decile_1_Simul.pdf')


% Figure with transitions of top decile
f=figure();
plot( fig_sample/1000, reshape(Tr_deciles_a_S(end,:,:),10,101))
yline(100*(Tr_deciles_a(end,:)), 'color', [0.8500 0.3250 0.0980])
ylim([0 ceil(max(max(reshape(Tr_deciles_a_S(end,:,:),10,101))))])
xlim([1 M_P.N_Panel/1000])
title('1st Decile Transitions','FontSize', 18)
xlabel('Sample Size: Thousands','FontSize',14)
f.Units = 'centimeters';
f.PaperUnits = 'centimeters';
f.PaperSize = f.Position(3:4);
f.PaperSize = f.Position(3:4)+0.2; % Add 0.1 cm of margin in each direction
saveas(gcf,'Fig_Folder/Tr_Decile_10_Simul.pdf')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 5 Year autocorrelation of consumption first quintile

av_c_q_S = zeros(2,fig_N) ; % Average Consumption in first quintile 
sd_c_q_S = zeros(2,fig_N) ; % Std Dev of Consumption in first quintile 
cor_c_q_S= zeros(fig_N,1) ; % Cor of Consumption in first quintile 



for i=1:fig_N
    % Fix current sample
    a_aux_0 = M_P.a_mat(1:fig_sample(i),8) ;
    c_aux_0 = M_P.c_mat(1:fig_sample(i),8) ;
    c_aux_T = M_P.c_mat(1:fig_sample(i),end) ;
    % Find index of first quintile 
    ind_q = find((0<=a_aux_0)&(a_aux_0<=deciles_a_S(3))) ;  %  Index of first quintile of assets in T-1
    % Compute moments
    av_c_q_S(1,i) = mean(c_aux_0(ind_q));
    sd_c_q_S(1,i) = std(c_aux_0(ind_q));
    av_c_q_S(2,i) = mean(c_aux_T(ind_q));
    sd_c_q_S(2,i) = std(c_aux_T(ind_q));
    cor_c_q_S(i) = corr(c_aux_0(ind_q), c_aux_T(ind_q) )   ;
end


% Figure
f=figure();
plot(fig_sample/1000, cor_c_q_S)
yline(mean(cor_c_q_S), 'color', [0.8500 0.3250 0.0980])
ylim([0 1])
xlim([1 M_P.N_Panel/1000])
title('Auto-Correlation of Consumption','FontSize', 18)
xlabel('Sample Size: Thousands','FontSize',14)
f.Units = 'centimeters';
f.PaperUnits = 'centimeters';
f.PaperSize = f.Position(3:4);
f.PaperSize = f.Position(3:4)+0.2; % Add 0.1 cm of margin in each direction
saveas(gcf,'Fig_Folder/Cons_Corr_Simul.pdf')



