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
med_eps = int64(round(Model.n_eps/2));
scatter(log(M_Aiyagari.eps_grid*par.w), 100*Gamma_eps, 200,'filled')
grid on
xlabel('(log) Labor Income','FontSize',18)
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
set(gca,'XTickLabelMode','auto')
text(log(M_Aiyagari.eps_grid(med_eps)*par.w)+0.15,2,sprintf('%3.1fk', round(M_Aiyagari.eps_grid(med_eps)*par.w,1)),'FontSize', 14, 'FontWeight','bold');
%text(log(M_Aiyagari.eps_grid(med_eps)*par.w)+0.3,2,sprintf('%3.1fk', round(M_Aiyagari.eps_grid(med_eps)*par.w,1)));
ylim([0 ceil(max(100*Gamma_eps/10))*10])
xlim([log(1) log(1000)])
xticks(log([1 5 10 50 100 500 1000]))
xticklabels({'1k', '5k', '10k', '50k', '100k', '500k', '1m' } )   
title('Labor Income Distribution','FontSize', 14)
hold on
xline(log(M_Aiyagari.eps_grid(med_eps)*par.w), '-k')
hold off
print('Fig_Folder/Distribution_Income_Shocks', '-dpdf', '-bestfit')
saveas(gcf,'Fig_Folder/Distribution_Income_Shocks.pdf')

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Asset Distribution

scatter(log(M_Aiyagari.a_grid_fine), 100*Gamma_a,'filled', 'blue')
grid on
xlabel('Log Assets','FontSize',18)
title('Asset Distribution','FontSize', 18)
ylim([0 ceil(max(100*Gamma_a))])
xlim([log(0.01) log(ceil(M_Aiyagari.a_grid(end)/1000)*1000)])
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
set(gca,'XTickLabelMode','auto')
xticks(log([1 10 100 1000 10000 100000]))
xticklabels({'1k', '10k', '100k', '1m', '10m', '100m' } )   
saveas(gcf,'Fig_Folder/Distribution_Wealth.pdf')

% Plot CDF

scatter(log(M_Aiyagari.a_grid_fine), 100*CDF_a, 'filled', 'blue')
grid on
xlabel('Log Assets','FontSize',18)
title('Cumulative Asset Distribution','FontSize', 18)
ylim([0 100])
xlim([log(0.1)  log(ceil(M_Aiyagari.a_grid(end)/1000)*1000)])
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
set(gca,'XTickLabelMode','auto')
xticks(log([1 10 100 1000 10000 100000]))
xticklabels({'1k', '10k', '100k', '1m', '10m', '100m' } )   
saveas(gcf,'Fig_Folder/Distribution_Wealth_CDF.pdf')


% Plot Pareto Tail (Above $1 Million)

ind = transpose(find(M_Aiyagari.a_grid_fine>=1000));
grid_1M = M_Aiyagari.a_grid_fine(ind) ;
Gamma_a_1M = Gamma_a(ind)/sum(Gamma_a(ind));
Gamma_a_1M = Gamma_a_1M/sum(Gamma_a_1M);
CCDF_1M = 1 - cumsum(Gamma_a_1M);
ind = transpose(find(CCDF_1M>= 1e-12));
P_coeff = dot(log(grid_1M(ind)/1000),log(grid_1M(ind)/1000)) \ dot(log(grid_1M(ind)/1000),log(CCDF_1M(ind)) ) ;

scatter(log(grid_1M(ind)/1000), log(CCDF_1M(ind)))
grid on
xlabel('Log Assets','FontSize',18)
title('Distribution Tail','FontSize', 18)
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
set(gca,'XTickLabelMode','auto')
ylim([floor(min(log(CCDF_1M(ind)))/4)*4  0])
xlim([log(1) log(ceil(M_Aiyagari.a_grid(end)/1000)*1)])
xticks(log([1 2 4 8 20 40 80]))
xticklabels({'1m', '2m', '4m', '8m', '20m', '40m' , '80m' } )   
text( -log(1.3)+mean(log(grid_1M(ind)/1000)) , mean(log(CCDF_1M(ind))) ,sprintf('α=%3.2f', round(P_coeff,2)),'FontSize', 14, 'FontWeight','bold');
hold on
str = '#D95319';
plot(log(grid_1M(1:end-1)/1000), P_coeff*log(grid_1M(1:end-1)/1000), '-', 'Color',  sscanf(str(2:end),'%2x%2x%2x',[1 3])/255 )
hold off
saveas(gcf,'Fig_Folder/Distribution_Wealth_Pareto.pdf')
clear str
% Plot Lorenz Curve

str = '#808080';
plot(1:100, 1:100, 'Color',  sscanf(str(2:end),'%2x%2x%2x',[1 3])/255)
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
set(gca,'XTickLabelMode','auto')
grid on
title('Lorenz Curve','FontSize', 18)
ylim([0 100])
xlim([0 100])
hold on
str = '#00688B';
plot(100*CDF_a, 100*Lorenz_a, 'Color',  sscanf(str(2:end),'%2x%2x%2x',[1 3])/255)
hold off
saveas(gcf,'Fig_Folder/Distribution_Wealth_Lorenz.pdf')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Age Profiles

% Labor Income
f = figure(); 
scatter(20:(19+M_Aiyagari.p.Max_Age), age_profile_y, 250, 'MarkerFaceColor', [0.3010 0.7450 0.9330]);
ylim([0 ceil(max(age_profile_y/10))*10])
title('Age Profile: Labor Income','FontSize', 38)
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
set(gca,'XTickLabelMode','auto')
grid on
xlabel('Age','FontSize',18)
ylabel('Thousands of Dollars','FontSize',18)
f.Units = 'centimeters';
f.PaperUnits = 'centimeters';
f.PaperSize = f.Position(3:4);
f.PaperSize = f.Position(3:4)+0.2; % Add 0.1 cm of margin in each direction
saveas(gcf,'Fig_Folder/Age_Profile_Y.pdf')


% Assets
f = figure(); 
scatter(20:(19+M_Aiyagari.p.Max_Age), age_profile_a, 250, 'MarkerFaceColor', [0.3010 0.7450 0.9330]);
ylim([0 ceil(max(age_profile_a/10))*10])
title('Age Profile: Assets','FontSize', 38)
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
set(gca,'XTickLabelMode','auto')
xticks(20:10:100)
grid on
xlabel('Age','FontSize',18)
ylabel('Thousands of Dollars','FontSize',18)
f.Units = 'centimeters';
f.PaperUnits = 'centimeters';
f.PaperSize = f.Position(3:4);
f.PaperSize = f.Position(3:4)+0.2; % Add 0.1 cm of margin in each direction
saveas(gcf,'Fig_Folder/Age_Profile_a.pdf')


% Savings Rate


f = figure(); 
scatter(20:(19+M_Aiyagari.p.Max_Age),  100*(age_profile_ap./age_profile_a-1) , 250, 'MarkerFaceColor', [0.3010 0.7450 0.9330]);
ylim([floor(min(min(100*(age_profile_ap./age_profile_a-1)/5)))*5   ...
    ceil(max(max(100*(age_profile_ap./age_profile_a-1)/5)))*5 ])
title('Age Profile: (agg) Savings Rate','FontSize', 38)
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
set(gca,'XTickLabelMode','auto')
xticks(20:10:100)
grid on
xlabel('Age','FontSize',18)
ylabel('Percentage Points','FontSize',18)
f.Units = 'centimeters';
f.PaperUnits = 'centimeters';
f.PaperSize = f.Position(3:4);
f.PaperSize = f.Position(3:4)+0.2; % Add 0.1 cm of margin in each direction
saveas(gcf,'Fig_Folder/Age_Profile_a.pdf')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Longitudinal Moments




% Get moments from histogram interation method

% Get moments from histogram method
run('CalculateMoments_Histogram.m')

% Get moments from Monte Carlo simulation
run("CalculateMoments_MonteCarlo.m")


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Solution Properties: Grids and Policy Functions



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot Grid and Fine Grid (Zoom at the bottom and the top)

f = figure(); 
subplot(2,2,[1,2]);
p1 = scatter(M_Aiyagari.a_grid, zeros(M_Aiyagari.n_a),'MarkerFaceColor', [0 0.4470 0.7410]);
xticks(0:2500:10000)
ylim([-0.25 1.25])
yticks([0 1])
yticklabels({'Coarse', 'Fine'})
xlabel('Assets (thousands of dollars)','FontSize',18)
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
set(gca,'XTickLabelMode','auto')
hold on
p2 = scatter(M_Aiyagari.a_grid, ones(M_Aiyagari.n_a),'MarkerFaceColor', [0.8500 0.3250 0.0980]);
hold off


subplot(2,2,3);
p3 = scatter(M_Aiyagari.a_grid, zeros(M_Aiyagari.n_a),'MarkerFaceColor', [0 0.4470 0.7410]);
xticks(0:2:10)
ylim([-0.25 1.25])
xlim([0 10])
yticks([0 1])
yticklabels({'Coarse', 'Fine'})
xlabel('Assets (thousands of dollars)','FontSize',18)
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
set(gca,'XTickLabelMode','auto')
hold on
p4 = scatter(M_Aiyagari.a_grid_fine, ones(M_Aiyagari.n_a_fine),'MarkerFaceColor', [0.8500 0.3250 0.0980]);
hold off


subplot(2,2,4);
p5 = scatter(M_Aiyagari.a_grid(1:end)./1000 , zeros(M_Aiyagari.n_a),'MarkerFaceColor', [0 0.4470 0.7410]);
xticks(0:2:10)
ylim([-0.25 1.25])
xlim([1 M_Aiyagari.a_max/1000])
yticks([0 1])
yticklabels({'Coarse', 'Fine'})
xlabel('Assets (millions of dollars)','FontSize',18)
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
set(gca,'XTickLabelMode','auto')
hold on
p6 = scatter(M_Aiyagari.a_grid_fine(1:end)./1000, ones(M_Aiyagari.n_a_fine),'MarkerFaceColor', [0.8500 0.3250 0.0980]);
hold off
f.Units = 'centimeters';
f.PaperUnits = 'centimeters';
f.PaperSize = f.Position(3:4);
f.PaperSize = f.Position(3:4)+0.2; % Add 0.1 cm of margin in each direction
saveas(gcf,'Fig_Folder/Asset_Grids_Layout.pdf')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot Saving Functions (median labor efficiency and interest rate)
f = figure();
plot(M_Aiyagari.a_grid_fine,M_Aiyagari.a_grid_fine, 'color', [.5 .5 .5])
title('Savings','FontSize', 38)
xlabel('Assets (thousands of dollars)','FontSize',18)
ylabel('Assets (thousands of dollars)','FontSize',18)
xlim([0 M_Aiyagari.a_max])
ylim([0 M_Aiyagari.a_max])
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
set(gca,'XTickLabelMode','auto')
hold on
plot(M_Aiyagari.a_grid_fine,M_Aiyagari.G_ap_fine(:,med_eps,ref_age),'color', [0 0.4470 0.7410])
hold off 
f.Units = 'centimeters';
f.PaperUnits = 'centimeters';
f.PaperSize = f.Position(3:4);
f.PaperSize = f.Position(3:4)+0.2; % Add 0.1 cm of margin in each direction
saveas(gcf,'Fig_Folder/Policy_Function_Savings_Level.pdf')



f = figure();
x = 100*(rdivide(M_Aiyagari.G_ap_fine(:,med_eps,ref_age),transpose(M_Aiyagari.a_grid_fine))-1) ;
plot( log(M_Aiyagari.a_grid_fine) ,   x ,'color', [0 0.4470 0.7410])
xticks(log([1 10 100 1000 10000 50000]))
xticklabels({'1k', '10k', '100k', '1m', '10m', '50m'} )   
xlim([log(0.8) log(M_Aiyagari.a_max)])
title('Savings Rate','FontSize', 38)
xlabel('(log) Assets','FontSize',18)
ylabel('Saving Rate (%)','FontSize',18)
yline(0,  'color', [.5 .5 .5])
ylim([-15 50])
f.Units = 'centimeters';
f.PaperUnits = 'centimeters';
f.PaperSize = f.Position(3:4);
f.PaperSize = f.Position(3:4)+0.2; % Add 0.1 cm of margin in each direction
saveas(gcf,'Fig_Folder/Policy_Function_Savings_Rate.pdf')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Euler Errors plots
Euler_perc_error = zeros(M_Aiyagari.n_a_fine,3);
 for i_a=1:M_Aiyagari.n_a_fine
        % Euler percentage error (lowest points in epsilon)
        Euler_perc_error(i_a,1) = Functions_ModelSolution.Euler_Error(1 ,ref_age,M_Aiyagari.a_grid_fine(i_a),M_Aiyagari);
        % Euler percentage error (median points in epsilon) 
        Euler_perc_error(i_a,2) = Functions_ModelSolution.Euler_Error(med_eps  ,ref_age,M_Aiyagari.a_grid_fine(i_a),M_Aiyagari);
        % Euler percentage error (higest points in epsilon)
        Euler_perc_error(i_a,3) = Functions_ModelSolution.Euler_Error(M_Aiyagari.n_eps,ref_age,M_Aiyagari.a_grid_fine(i_a),M_Aiyagari);
 end

 f=figure();
 subplot(3,1,1)
p1_euler = scatter(log(M_Aiyagari.a_grid_fine), Euler_perc_error(:,1));
ylabel('Percentage','FontSize',12)
xlim([log(0.8) log(M_Aiyagari.a_max + 0.2)]);
xticks(log([1 10 100 1000 10000 50000]))
xticklabels({'1k', '10k', '100k', '1m', '10m', '50m'} )   
ylim([min(Euler_perc_error(:,1))  max(Euler_perc_error(:,1)) ])

 subplot(3,1,2)
p2_euler = scatter(log(M_Aiyagari.a_grid_fine), Euler_perc_error(:,2));
ylabel('Percentage','FontSize',12)
xlim([log(0.8) log(M_Aiyagari.a_max + 0.2)]);
xticks(log([1 10 100 1000 10000 50000]))
xticklabels({'1k', '10k', '100k', '1m', '10m', '50m'} )   
ylim([min(Euler_perc_error(:,2))  max(Euler_perc_error(:,2)) ])

 subplot(3,1,3)
 p3_euler = scatter(log(M_Aiyagari.a_grid_fine), Euler_perc_error(:,3));
ylabel('Percentage','FontSize',12)
xlim([log(0.8) log(M_Aiyagari.a_max + 0.2)]);
xticks(log([1 10 100 1000 10000 50000]))
xticklabels({'1k', '10k', '100k', '1m', '10m', '50m'} )   
ylim([min(Euler_perc_error(:,3))  max(Euler_perc_error(:,3)) ])

f.Units = 'centimeters';
f.PaperUnits = 'centimeters';
f.PaperSize = f.Position(3:4);
f.PaperSize = f.Position(3:4)+0.2; % Add 0.1 cm of margin in each direction
saveas(gcf,'Fig_Folder/Euler_Error_Layout.pdf')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
