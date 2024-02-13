%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Graphs and Tables for Draft
%%% Loads results generated by Draft_Results

color_vec_H = {'k', [0.53, 0.81, 0.98], [0.39, 0.58, 0.93], [0.25, 0.41, 0.88]};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load results from csv files

H_grid_size = [250, 250, 500, 1000];
n_H = length(H_grid_size);
pct_list = [90; 95; 99; 99.9; 99.99];
S_sample = [50000, 100000, 250000, 500000];
N_S = length(S_sample);
age_0_Wealth_Profile = 26;
age_0_Wealth_Corr = 16;
age_T_Wealth_Corr = 36;

H_G_timed = dlmread(fullfile(Hist_Folder, 'H_G_timed.csv'), ',', 0, 0);
H_G_bytes = dlmread(fullfile(Hist_Folder, 'H_G_bytes.csv'), ',', 0, 0);
H_M_timed = reshape(dlmread(fullfile(Hist_Folder, 'H_M_timed.csv'), ',', 0, 0), n_H, 3);
H_M_bytes = reshape(dlmread(fullfile(Hist_Folder, 'H_M_bytes.csv'), ',', 0, 0), n_H, 3);

H_Wealth_Profile_NB = reshape(dlmread(fullfile(Hist_Folder, 'H_Wealth_Profile_NB.csv'), ',', 0, 0), p.Max_Age, 6, N_S);
H_Wealth_Profile_45 = reshape(dlmread(fullfile(Hist_Folder, 'H_Wealth_Profile_45.csv'), ',', 0, 0), p.Max_Age - age_0_Wealth_Profile + 1, 6, N_S);
H_Wealth_Corr = reshape(dlmread(fullfile(Hist_Folder, 'H_Wealth_Corr.csv'), ',', 0, 0), n_H, 2);

S_M_timed = reshape(dlmread(fullfile(MC_Folder, 'S_M_timed.csv'), ',', 0, 0), N_S, 6);
S_M_bytes = reshape(dlmread(fullfile(MC_Folder, 'S_M_bytes.csv'), ',', 0, 0), N_S, 4);

S_Wealth_Profile_NB = reshape(dlmread(fullfile(MC_Folder, 'S_Wealth_Profile_NB.csv'), ',', 0, 0), p.Max_Age, 6, N_S);
S_Wealth_Profile_45 = reshape(dlmread(fullfile(MC_Folder, 'S_Wealth_Profile_45.csv'), ',', 0, 0), p.Max_Age, 6, N_S);
S_Wealth_Corr = reshape(dlmread(fullfile(MC_Folder, 'S_Wealth_Corr.csv'), ',', 0, 0), N_S, 2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Wealth Profiles Figures 


y_tick_pr = [0, 250, 500, 750, 1000, 1250, 1500, 1750];
y_label_pr = {'$0m', '$0.25m', '$0.50m', '$0.75m', '$1.00m', '$1.25m', '$1.50m', '$1.75m'};
S_label = {'100k', '250k', '350k', '500k'};

for i = 1:N_S
    figure;
    set(gcf, 'Visible', 'off');  % To prevent the figures from being displayed

    scatter(0, 0, 'Marker', 'd', 'MarkerSize', 15, 'MarkerEdgeColor', 'white', 'DisplayName', 'Ave. Wealth');
    hold on;
    scatter(0, 0, 'o', 'MarkerSize', 15, 'MarkerEdgeColor', 'white', 'DisplayName', 'p99 Wealth');
    scatter(0, 0, 'd', 'MarkerSize', 15, 'MarkerEdgeColor', 'white', 'DisplayName', 'p99.9 Wealth');
    scatter(0, 0, 'Marker', 's', 'MarkerSize', 15, 'MarkerEdgeColor', 'cornflowerblue', 'DisplayName', 'Histogram, N=500');
    scatter(0, 0, 'Marker', 's', 'MarkerSize', 15, 'MarkerEdgeColor', 'orangered2', 'DisplayName', ['Monte-Carlo, N=' S_label{i}]);

    scatter(45:(19 + p.Max_Age), H_Wealth_Profile_NB(age_0_Wealth_Profile:end, end, 2), 'Marker', 'd', 'MarkerSize', 4, 'MarkerFaceColor', 'cornflowerblue', 'MarkerEdgeColor', 'cornflowerblue', 'DisplayName', '');
    scatter(45:(19 + p.Max_Age), S_Wealth_Profile_NB(age_0_Wealth_Profile:end, end, i), 'Marker', 'd', 'MarkerSize', 4, 'MarkerFaceColor', 'orangered2', 'MarkerEdgeColor', 'grey25', 'DisplayName', '');

    scatter(45:(19 + p.Max_Age), H_Wealth_Profile_NB(age_0_Wealth_Profile:end, 3, 2), 'o', 'MarkerSize', 4.25, 'MarkerFaceColor', 'cornflowerblue', 'MarkerEdgeColor', 'cornflowerblue', 'DisplayName', '');
    scatter(45:(19 + p.Max_Age), S_Wealth_Profile_NB(age_0_Wealth_Profile:end, 3, i), 'o', 'MarkerSize', 4.25, 'MarkerFaceColor', 'orangered2', 'MarkerEdgeColor', 'grey25', 'DisplayName', '');

    scatter(45:(19 + p.Max_Age), H_Wealth_Profile_NB(age_0_Wealth_Profile:end, 4, 2), 'd', 'MarkerSize', 4.25, 'MarkerFaceColor', 'cornflowerblue', 'MarkerEdgeColor', 'cornflowerblue', 'DisplayName', '');
    scatter(45:(19 + p.Max_Age), S_Wealth_Profile_NB(age_0_Wealth_Profile:end, 4, i), 'd', 'MarkerSize', 4.25, 'MarkerFaceColor', 'orangered2', 'MarkerEdgeColor', 'grey25', 'DisplayName', '');

    ylim([0, ceil(max(max(S_Wealth_Profile_NB(:, 4, 1) / 250)) * 250)]);
    xlabel('Age', 'FontSize', 18);
    ylabel('Thousands of Dollars', 'FontSize', 18);
    xticks(45:10:100);
    xlim([45, 100]);

    legend('Location', 'northeastoutside', 'FontSize', 9);
    saveas(gcf, fullfile(Fig_Folder, ['Draft_Wealth_Profile_45_' num2str(S_sample(i) / 1000) 'k.pdf']));

    close gcf;  % Close the figure to avoid displaying it
end

% Comparing Histograms
alpha_vec = linspace(0.85, 0.95, 4);
ms_vec = linspace(4.50, 3.25, 4);

figure;
set(gcf, 'Visible', 'off');

scatter(0, 0, 'Marker', 'd', 'MarkerSize', 15, 'MarkerEdgeColor', 'white', 'DisplayName', 'Ave. Wealth');
hold on;
scatter(0, 0, 'o', 'MarkerSize', 15, 'MarkerEdgeColor', 'white', 'DisplayName', 'p99 Wealth');
scatter(0, 0, 'd', 'MarkerSize', 15, 'MarkerEdgeColor', 'white', 'DisplayName', 'p99.9 Wealth');

for i = [2, 3, 4]
    scatter(0, 0, 'Marker', 's', 'MarkerSize', 15, 'MarkerEdgeColor', color_vec_H{i}, 'DisplayName', ['Histogram, N=' num2str(H_grid_size(i))]);

    scatter(45:(19 + p.Max_Age), H_Wealth_Profile_NB(age_0_Wealth_Profile:end, end, i), 'Marker', 'd', 'MarkerSize', ms_vec(i), 'MarkerFaceAlpha', alpha_vec(i), 'MarkerEdgeColor', color_vec_H{i}, 'DisplayName', '');
    scatter(45:(19 + p.Max_Age), H_Wealth_Profile_NB(age_0_Wealth_Profile:end, 3, i), 'o', 'MarkerSize', ms_vec(i), 'MarkerFaceAlpha', alpha_vec(i), 'MarkerEdgeColor', color_vec_H{i}, 'DisplayName', '');
    scatter(45:(19 + p.Max_Age), H_Wealth_Profile_NB(age_0_Wealth_Profile:end, 4, 2), 'd', 'MarkerSize', ms_vec(i), 'MarkerFaceAlpha', alpha_vec(i), 'MarkerEdgeColor', color_vec_H{i}, 'DisplayName', '');
end

ylim([0, ceil(max(max(S_Wealth_Profile_NB(:, 4, 1) / 250)) * 250]);
yticks(y_tick_pr);
xlabel('Age', 'FontSize', 18);
ylabel('Thousands of Dollars', 'FontSize', 18);
xticks(45:10:100);
xlim([45, 100]);

legend('Location', 'northeastoutside', 'FontSize', 9);

saveas(gcf, fullfile(Fig_Folder, 'Draft_Wealth_Profile_45_Hist.pdf'));
close gcf;

% Figure for legend
step_1 = 2.6;
step_2 = 3.5;
x_vec_m = [(0.50:step_1:(0.5 + step_1 * 3))', (0.5 + step_1 * 4):step_2:(0.5 + step_1 * 4 + step_2 * 2)'];
x_vec_a = x_vec_m + 0.30;

figure;
set(gca, 'Visible', 'off');

scatter(x_vec_m(1), 1.0, 'Marker', 'd', 'MarkerSize', 4.5, 'MarkerEdgeColor', 'white');
annotation('textbox', [x_vec_a(1), 1.0, 0, 0], 'String', 'Ave. Wealth', 'HorizontalAlignment', 'left', 'FontSize', 5);

for j = 2:7
    scatter(x_vec_m(j), 1.0, 'Marker', 's', 'MarkerSize', 4.5, 'MarkerEdgeColor', color_vec_H{j - 1});
    annotation('textbox', [x_vec_a(j), 1.0, 0, 0], 'String', ['Histogram, N=' num2str(H_grid_size(j - 1))], 'HorizontalAlignment', 'left', 'FontSize', 5);
end

annotation('textbox', [x_vec_m(end) + step_2 + 0.10, 0.97, 0, 0], 'String', 'Monte-Carlo', 'HorizontalAlignment', 'left', 'FontSize', 5);
annotation('line', [0, 0, x_vec_m(end) + step_2 + 0.10, x_vec_m(end) + step_2 + 0.10, 0], [0.97, 1.03, 1.03, 0.97, 0.97], 'Color', [0.75, 0.75, 0.75], 'LineWidth', 0.5);

xlim([-0.5, x_vec_m(end) + step_2 + 1.0]);
ylim([0.5, 1.5]);

saveas(gcf, fullfile(Fig_Folder, 'Draft_Wealth_Profile_45_Legend.pdf'));
close gcf;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Autocorrelation

Mat_Corr = [
    "Wealth Autocorrelation 35-65", "", "", "", "";
    "Histogram", "", "", "", "";
    "Grid Size", H_grid_size;
    "Ave. a(35)", H_Wealth_Profile_NB(age_0_Wealth_Corr, end, :)';
    "p99  a(35)", H_Wealth_Profile_NB(age_0_Wealth_Corr, 3, :)';
    "Ave. a(65)", H_Wealth_Profile_NB(age_T_Wealth_Corr, end, :)';
    "p99  a(65)", H_Wealth_Profile_NB(age_T_Wealth_Corr, 3, :)';
    "Wealth Corr 35-40", 100 * H_Wealth_Corr(:, 1)';
    "Wealth Corr 35-50", 100 * H_Wealth_Corr(:, 2)';
    "-" "-" "-" "-" "-";
    "Model Time", H_Γ_timed';
    "Profile Time", H_M_timed(:, 1)';
    "Corr Time 35-45", H_M_timed(:, 2)';
    "Corr Time 35-45", H_M_timed(:, 3)';
    "-" "-" "-" "-" "-";
    "-" "-" "-" "-" "-";
    "-" "-" "-" "-" "-";
    "Simulation", "" "" "" "";
    "Sample Size (k)", int32(S_sample / 1000)';
    "Ave. a(35)", S_Wealth_Profile_NB(age_0_Wealth_Corr, end, :)';
    "p99  a(35)", S_Wealth_Profile_NB(age_0_Wealth_Corr, 3, :)';
    "Ave. a(65)", S_Wealth_Profile_NB(age_T_Wealth_Corr, end, :)';
    "p99  a(65)", S_Wealth_Profile_NB(age_T_Wealth_Corr, 3, :)';
    "Wealth Corr 35-40", 100 * S_Wealth_Corr(:, 1)';
    "Wealth Corr 35-50", 100 * S_Wealth_Corr(:, 2)';
    "-" "-" "-" "-" "-";
    "Simulation Time", S_M_timed(:, 1)';
    "Cohort Simul Time", S_M_timed(:, 6)';
    "Profile Time", S_M_timed(:, 2)';
    "Corr Time 35-45", S_M_timed(:, 3)';
    "Corr Time 35-55", S_M_timed(:, 4)';
    "-" "-" "-" "-" "-";
];

% Write to CSV file
csvwrite(fullfile(Fig_Folder, 'Table_Auto_Corr.csv'), Mat_Corr);

