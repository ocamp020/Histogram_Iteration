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

% Plot Income Distribution

med_eps = int64(round(Model.n_eps/2));
f=figure();
scatter(log(M_Aiyagari.eps_grid*par.w), 100*Gamma_eps, 200,'filled')
grid on
xlabel('(log) Labor Income','FontSize',18)
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
f.Units = 'centimeters';
f.PaperUnits = 'centimeters';
f.PaperSize = f.Position(3:4);
f.PaperSize = f.Position(3:4)+0.2; % Add 0.1 cm of margin in each direction
saveas(gcf,'Fig_Folder/Distribution_Income.pdf')



% Plot Return Distribution
med_eps = int64(round(Model.n_eps/2));
f=figure();
scatter(M_Aiyagari.zeta_grid*par.r*100, 100*Gamma_zeta, 200,'filled')
text(M_Aiyagari.zeta_grid*par.r*100,100*Gamma_zeta, cellstr(num2str(transpose(round(M_Aiyagari.zeta_grid*par.r*100,2)))), ...
    'VerticalAlignment','top','HorizontalAlignment','left')
grid on
xlabel('Percentage Points','FontSize',18)
title('Return Distribution','FontSize', 14)
xlim([0 ceil(max(100*M_Aiyagari.zeta_grid*p.r/10))*10])
ylim([0 ceil(max(100*Gamma_zeta/10))*10])
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
set(gca,'XTickLabelMode','auto')
f.Units = 'centimeters';
f.PaperUnits = 'centimeters';
f.PaperSize = f.Position(3:4);
f.PaperSize = f.Position(3:4)+0.2; % Add 0.1 cm of margin in each direction
saveas(gcf,'Fig_Folder/Distribution_Return.pdf')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot asset distribution
f=figure();
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
f.Units = 'centimeters';
f.PaperUnits = 'centimeters';
f.PaperSize = f.Position(3:4);
f.PaperSize = f.Position(3:4)+0.2; % Add 0.1 cm of margin in each direction
saveas(gcf,'Fig_Folder/Distribution_Wealth.pdf')


% Plot CDF
f=figure();
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
f.Units = 'centimeters';
f.PaperUnits = 'centimeters';
f.PaperSize = f.Position(3:4);
f.PaperSize = f.Position(3:4)+0.2; % Add 0.1 cm of margin in each direction
saveas(gcf,'Fig_Folder/Distribution_Wealth_CDF.pdf')




% Plot Pareto Tail (Above $5 Million)

ind = transpose(find(M_Aiyagari.a_grid_fine>=5000));
grid_1M = M_Aiyagari.a_grid_fine(ind) ;
Gamma_a_1M = Gamma_a(ind)/sum(Gamma_a(ind));
Gamma_a_1M = Gamma_a_1M/sum(Gamma_a_1M);
CCDF_1M = 1 - cumsum(Gamma_a_1M);
ind = transpose(find(CCDF_1M>= 1e-12));
P_coeff = dot(log(grid_1M(ind)/5000),log(grid_1M(ind)/5000)) \ dot(log(grid_1M(ind)/5000),log(CCDF_1M(ind)) ) ;
scatter(log(grid_1M(ind)/5000), log(CCDF_1M(ind)))
grid on
xlabel('Log Assets','FontSize',18)
title('Distribution Tail','FontSize', 18)
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
set(gca,'XTickLabelMode','auto')
ylim([floor(min(log(CCDF_1M(ind)))/4)*4  0])
xlim([log(1) log(ceil(M_Aiyagari.a_grid(end)/5000)*1)])
xticks(log([1 2 4 8 16 40 80]))
xticklabels({'5m', '10m', '20m', '40m', '80m' } )   
text( -log(1.3)+mean(log(grid_1M(ind)/5000)) , mean(log(CCDF_1M(ind))) ,sprintf('α=%3.2f', round(P_coeff,2)),'FontSize', 14, 'FontWeight','bold');
hold on
str = '#D95319';
plot(log(grid_1M(1:end-1)/5000), P_coeff*log(grid_1M(1:end-1)/5000), '-', 'Color',  sscanf(str(2:end),'%2x%2x%2x',[1 3])/255 )
hold off
saveas(gcf,'Fig_Folder/Distribution_Wealth_Pareto.pdf')
clear str

% Plot Lorenz Curve

str = '#808080';
plot(1:100, 1:100, 'Color',  sscanf(str(2:end),'%2x%2x%2x',[1 3])/255,'LineWidth',2.0)
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
set(gca,'XTickLabelMode','auto')
grid on
title('Lorenz Curve','FontSize', 18)
ylim([0 100])
xlim([0 100])
hold on
str = '#00688B';
plot(100*CDF_a, 100*Lorenz_a, 'Color',  sscanf(str(2:end),'%2x%2x%2x',[1 3])/255, 'LineWidth',2.0)
hold off
saveas(gcf,'Fig_Folder/Distribution_Wealth_Lorenz.pdf')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



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


f = figure(); 
subplot(2,2,[1,2]);
p1 = scatter(M_Aiyagari.a_grid, zeros(M_Aiyagari.n_a),'MarkerFaceColor', [0 0.4470 0.7410]);
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
% Plot Saving Functions (median labor efficiency and interest rate)

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
plot(M_Aiyagari.a_grid_fine,M_Aiyagari.G_ap_fine(:,med_eps,med_zeta),'color', [0 0.4470 0.7410])
hold off 
f.Units = 'centimeters';
f.PaperUnits = 'centimeters';
f.PaperSize = f.Position(3:4);
f.PaperSize = f.Position(3:4)+0.2; % Add 0.1 cm of margin in each direction
saveas(gcf,'Fig_Folder/Policy_Function_Savings_Level.pdf')



f = figure();
x = 100*(rdivide(M_Aiyagari.G_ap_fine(:,med_eps,med_zeta),transpose(M_Aiyagari.a_grid_fine))-1) ;
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
% Euler Errors plots
Euler_perc_error = zeros(M_Aiyagari.n_a_fine,3);
 for i_a=1:M_Aiyagari.n_a_fine
        % Euler percentage error (lowest points in epsilon)
        Euler_perc_error(i_a,1) = Functions_ModelSolution.Euler_Error(1 ,1,M_Aiyagari.a_grid_fine(i_a),M_Aiyagari);
        % Euler percentage error (median points in epsilon) 
        Euler_perc_error(i_a,2) = Functions_ModelSolution.Euler_Error(med_eps  ,med_zeta,M_Aiyagari.a_grid_fine(i_a),M_Aiyagari);
        % Euler percentage error (higest points in epsilon)
        Euler_perc_error(i_a,3) = Functions_ModelSolution.Euler_Error(M_Aiyagari.n_eps,M_Aiyagari.n_zeta,M_Aiyagari.a_grid_fine(i_a),M_Aiyagari);
 end




f=figure();
% Plot lowest epsilon, lowest zeta
subplot(3,1,1)
p1_euler = scatter(log(M_Aiyagari.a_grid_fine), Euler_perc_error(:,1));
ylabel('Percentage','FontSize',12)
xlim([log(0.8) log(M_Aiyagari.a_max + 0.2)]);
xticks(log([1 10 100 1000 10000 50000]))
xticklabels({'1k', '10k', '100k', '1m', '10m', '50m'} )   
ylim([min(Euler_perc_error(:,1))  max(Euler_perc_error(:,1)) ])

% Plot median epsilon, median zeta
 subplot(3,1,2)
p2_euler = scatter(log(M_Aiyagari.a_grid_fine), Euler_perc_error(:,2));
ylabel('Percentage','FontSize',12)
xlim([log(0.8) log(M_Aiyagari.a_max + 0.2)]);
xticks(log([1 10 100 1000 10000 50000]))
xticklabels({'1k', '10k', '100k', '1m', '10m', '50m'} )   
ylim([min(Euler_perc_error(:,2))  max(Euler_perc_error(:,2)) ])

% Plot highest epsilon, highest zeta
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



