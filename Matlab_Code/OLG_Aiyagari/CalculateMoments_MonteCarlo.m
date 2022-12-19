%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute moments using the Monte-Carlo method 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run Simulations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate structure for the panel 

%M_P = struct('N_Panel', 500000 , 'T_Panel', 10, 'T_Simul', 1500, 'N_Min', 1000, 'Simul_tol', 1e-2, 'rng_seed', 297835398); 
% Panel Output
%M_P.a_mat    = zeros([M_P.N_Panel, M_P.T_Panel]);
%M_P.c_mat    =  zeros([M_P.N_Panel, M_P.T_Panel]);
%M_P.eps_mat  = zeros([M_P.N_Panel, M_P.T_Panel]);
%M_P.h_mat    =  zeros([M_P.N_Panel, M_P.T_Panel]);
%M_P.t_mat    = zeros([M_P.N_Panel,1]);

M_P = Functions_MonteCarlo.Model_Panel(500000, 10, 1500, 1000, 1e-2, 297835398);

M_P = Functions_MonteCarlo.Simulate_Panel(M_Aiyagari, M_P) ;

fig_sample = [1000 ; transpose(5000:5000:M_P.N_Panel)];
fig_N      = length(fig_sample)         ;

% Save results

writematrix(M_P.a_mat,'File_Folder/S_a_mat.csv')
writematrix(M_P.c_mat,'File_Folder/S_c_mat.csv')
writematrix(M_P.eps_mat,'File_Folder/S_e_mat.csv')
writematrix(M_P.h_mat,'File_Folder/S_h_mat.csv')



% Model Cohort
M_C_45 = Functions_MonteCarlo.Model_Cohort(500000,  M_Aiyagari.p.Max_Age-25,297835398 );
%M_C_45 =  struct('N_Panel', 500000, 'T_Panel', M_Aiyagari.p.Max_Age-25,  'rng_seed', 297835398) ;
%M_C_45.a_mat    = zeros([M_C_45.N_Panel, M_C_45.T_Panel]);
%M_C_45.c_mat    =  zeros([M_C_45.N_Panel, M_C_45.T_Panel]);
%M_C_45.eps_mat  = zeros([M_C_45.N_Panel, M_C_45.T_Panel]);
%M_C_45.h_mat    =  zeros([M_C_45.N_Panel, M_C_45.T_Panel]);
%M_C_45.t_mat    = zeros([M_C_45.N_Panel,1]);


M_C_45 = Functions_MonteCarlo.Simulate_Cohort(M_Aiyagari, M_C_45, 26);


M_C_35 = Functions_MonteCarlo.Model_Cohort(500000,  M_Aiyagari.p.Max_Age-15,297835398 );
%M_C_35 =  struct('N_Panel', 500000, 'T_Panel', M_Aiyagari.p.Max_Age-15,  'rng_seed', 297835398) ;
%M_C_35.a_mat    = zeros([M_C_35.N_Panel, M_C_35.T_Panel]);
%M_C_35.c_mat    =  zeros([M_C_35.N_Panel, M_C_35.T_Panel]);
%M_C_35.eps_mat  = zeros([M_C_35.N_Panel, M_C_35.T_Panel]);
%M_C_35.h_mat    =  zeros([M_C_35.N_Panel, M_C_35.T_Panel]);
%M_C_35.t_mat    = zeros([M_C_35.N_Panel,1]);

M_C_35 = Functions_MonteCarlo.Simulate_Cohort(M_Aiyagari, M_C_35, 16);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cohort evolution of average and 90-10 percentiles of wealth 

disp("---------------------------------------------------------")
disp("Average and 90-10 Percentiles for Cohorts")
% Newborns
sample_vec = [10000, 50000, 100000, 500000]; 
pct_list = [10;50;90;99];
av_a_S  = zeros(M_Aiyagari.p.Max_Age,4);
pct_S  = zeros(4,M_Aiyagari.p.Max_Age,4); 
for i=1:4
    a_sample = M_P.a_mat(1:sample_vec(i),end);
    h_sample = M_P.h_mat(1:sample_vec(i),end);
    for age=1:M_Aiyagari.p.Max_Age
        if sum(h_sample==age)>0
            av_a_S(age,i) = mean(a_sample(h_sample==age));
            pct_S(:,age,i) = prctile(a_sample(h_sample==age), pct_list);
        else
            fprintf('Error in sample at i=%i, sample=%i and age=%i: (%i) Observations \n', i, sample_vec(i), age,sum(h_sample==age))
            av_a_S(age,i) = nan;
            pct_S(:,age,i) = nan;
        end
    end
end

% Figure Average
f=figure();
scatter( 20:(19+M_Aiyagari.p.Max_Age) , age_profile_a,'MarkerFaceColor', [0 0.4470 0.7410]);
ylim([0 ceil(max(age_profile_a/10))*10])
title('Age Profile: Average Assets','FontSize', 18)
xticks(20:10:100)
ylabel('Thousands of dollars','FontSize',14)
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
set(gca,'XTickLabelMode','auto')
hold on
scatter( 20:(19+M_Aiyagari.p.Max_Age), av_a_S(:,1),'MarkerFaceColor', [0.8500 0.3250 0.0980]);
scatter( 20:(19+M_Aiyagari.p.Max_Age), av_a_S(:,3),'MarkerFaceColor', [1 0 0]);
scatter( 20:(19+M_Aiyagari.p.Max_Age), av_a_S(:,4),'MarkerFaceColor', [0 1 0]);
hold off
f.Units = 'centimeters';
f.PaperUnits = 'centimeters';
f.PaperSize = f.Position(3:4);
f.PaperSize = f.Position(3:4)+0.2; % Add 0.1 cm of margin in each direction
saveas(gcf,'Fig_Folder/Age_Profile_a_Simul.pdf')




% Figure Percentiles
f=figure();
scatter( 20:(19+M_Aiyagari.p.Max_Age) , pct_C_age_1(1,:) ,'MarkerFaceColor', [0 0.4470 0.7410]);
ylim([0 ceil(max(pct_C_age_1(2,:)/10))*10])
title('Age Profile: 10th and 90th Percentiles of Assets','FontSize', 18)
ylabel('Thousands of dollars','FontSize',14)
xlabel('Age','FontSize',14)
xticks(20:10:100)
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
set(gca,'XTickLabelMode','auto')
hold on
scatter( 20:(19+M_Aiyagari.p.Max_Age) , pct_S(1,:,1) );

scatter( 20:(19+M_Aiyagari.p.Max_Age) , pct_S(1,:,3));
scatter( 20:(19+M_Aiyagari.p.Max_Age) , pct_S(1,:,4));

% 90 Percentile
scatter( 20:(19+M_Aiyagari.p.Max_Age) , pct_C_age_1(2,:),'diamond' ,'MarkerFaceColor', [0 0.4470 0.7410]);
scatter( 20:(19+M_Aiyagari.p.Max_Age) , pct_S(3,:,1), 'diamond' );

scatter( 20:(19+M_Aiyagari.p.Max_Age) , pct_S(3,:,3), 'diamond' );
scatter( 20:(19+M_Aiyagari.p.Max_Age) , pct_S(3,:,4), 'diamond' );
hold off
f.Units = 'centimeters';
f.PaperUnits = 'centimeters';
f.PaperSize = f.Position(3:4);
f.PaperSize = f.Position(3:4)+0.2; % Add 0.1 cm of margin in each direction
saveas(gcf,'Fig_Folder/Age_Profile_a_pct_Simul.pdf')




% 45-year olds with more than median income
age_0=26;
sample_vec = [10000, 50000, 100000, 500000]  ; 
pct_list   = [10;50;90;99]                   ;
pct_2_S    = zeros(4,M_Aiyagari.p.Max_Age-(age_0-1),4) ; 
for i=1:4
    a_sample    = M_C_45.a_mat(1:sample_vec(i),:)    ;
    h_sample    = M_C_45.h_mat(1:sample_vec(i),:)    ;
    eps_sample    = M_C_45.eps_mat(1:sample_vec(i),:)    ;
    ind_45      = find(eps_sample(:,1)>=med_eps) ;
    a_sample    = a_sample(ind_45,:);
    h_sample    = h_sample(ind_45,:);  
    for age=1:(M_Aiyagari.p.Max_Age-(age_0-1))
        ind_alive = find(h_sample(:,age)==1);
        if sum(ind_alive)>0
            pct_2_S(:,age,i) = prctile(a_sample(ind_alive,age),pct_list);
        else
            pct_2_S(:,age,i) = nan;
        end
    end
end


% Figure Percentiles for 45 year old 
f=figure();
scatter( 45:(19+M_Aiyagari.p.Max_Age) , pct_C_age_2(2,:) ,'MarkerFaceColor', [0 0.4470 0.7410]);
ylim([0 ceil(max(pct_C_age_2(2,:)/10))*10])
title('Age Profile: 10th and 90th Percentiles of Assets','FontSize', 18)
ylabel('Thousands of dollars','FontSize',14)
xlabel('Age','FontSize',14)
xticks(45:10:100)
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
set(gca,'XTickLabelMode','auto')
hold on
scatter( 45:(19+M_Aiyagari.p.Max_Age) , pct_2_S(1,:,1) );

scatter( 45:(19+M_Aiyagari.p.Max_Age) , pct_2_S(1,:,3));
scatter( 45:(19+M_Aiyagari.p.Max_Age) , pct_2_S(1,:,4));

% 90 Percentile
scatter( 45:(19+M_Aiyagari.p.Max_Age) , pct_C_age_2(2,:),'diamond' ,'MarkerFaceColor', [0 0.4470 0.7410]);
scatter( 45:(19+M_Aiyagari.p.Max_Age) , pct_2_S(3,:,1), 'diamond' );

scatter( 45:(19+M_Aiyagari.p.Max_Age) , pct_2_S(3,:,3), 'diamond' );
scatter( 45:(19+M_Aiyagari.p.Max_Age) , pct_2_S(3,:,4), 'diamond' );
hold off
f.Units = 'centimeters';
f.PaperUnits = 'centimeters';
f.PaperSize = f.Position(3:4);
f.PaperSize = f.Position(3:4)+0.2; % Add 0.1 cm of margin in each direction
saveas(gcf,'Fig_Folder/Age_Profile_a_pct_45_Simul.pdf')





   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auto-correlation of wealth Ages 35 and 65 (conditional on survival)

disp("\n===============================================")
disp("Autocorr of Wealth: ages 35-65")

age_0 = 16; %Start agents at age 35
n_H          = 65-(35-1); % Simulate for n_H years 
cor_a_3565_S = zeros(fig_N,1); % Cor of Consumption in first quintile 
for i=1:fig_N
    % Fix current sample
    a_aux_0 = M_C_35.a_mat(1:fig_sample(i),1);
    a_aux_T = M_C_35.a_mat(1:fig_sample(i),n_H);
    % Compute moments
    cor_a_3565_S(i) = corr(a_aux_0, a_aux_T) ;
end


% Figure
f=figure();
plot( fig_sample/1000, cor_a_3565_S,'LineWidth',2 )
yline(mean(cor_a_3565_S), 'color', [0.8500 0.3250 0.0980],'LineWidth',1)
ylim([0 1])
xlim([1 M_P.N_Panel/1000])
title('Auto-Correlation of Wealth, ages 35-65','FontSize', 18)
xlabel('Sample Size: Thousands','FontSize',14)
f.Units = 'centimeters';
f.PaperUnits = 'centimeters';
f.PaperSize = f.Position(3:4);
f.PaperSize = f.Position(3:4)+0.2; % Add 0.1 cm of margin in each direction
saveas(gcf,'Fig_Folder/A_Corr_3565_Simul.pdf')





disp("===============================================\n")

