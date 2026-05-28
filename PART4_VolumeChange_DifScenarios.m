%% Figure S3H: Volume loss for different membrane elasticities
% This panel quantifies how the maximum volume loss depends on the
% membrane elastic modulus \tilde{k}.

k_all = [0.1 0.25 0.5];

for kk = 1:length(k_all)
    CaseName = strcat('Vloss_Elasticity_', num2str(k_all(kk)));
    load(strcat('Data_for_plot\', CaseName, '.mat'));
    
    % Fractional volume change relative to the initial state
    V_all = Geometry{1};
    dV_all_perc = (V_all - V_all(1,:)) ./ V_all(1,:);
    
    % Maximum volume loss across all cells and time points
    Tension_Vloss_all_max(kk) = min(dV_all_perc(:));
end

% Plot as bar chart
figure;

% Use the default MATLAB axis color
co = get(groot, 'defaultAxesColorOrder');
c1 = co(1,:);

b = bar(Tension_Vloss_all_max(end:-1:1)); hold on;

% Bar appearance
b.FaceColor = c1;
b.EdgeColor = c1;
b.LineWidth = 1.5;

% Axis formatting
set(gca, 'FontSize', 18, 'FontName', 'Arial');
xticks(1:length(k_all));
xticklabels({'0.5','0.25','0.1'});

xlabel('$\tilde{k}$', 'Interpreter', 'latex');
ylabel('Volume loss (%)');

axis square;
box on;

%% Figure 5E: Volume loss versus pressure differential for different mean initial curvatures
% This panel examines how the response amplitude varies with the initial
% hydraulic pressure difference across the perturbed cells for simulations
% with different mean initial curvatures.

clc; clear; close all;

% File settings
dataFolder = 'Data_for_plot\Global_Varying_mean_initial_curvature\';
fileList = dir(fullfile(dataFolder, 'C0_*.mat'));
nCase = numel(fileList);

if nCase == 0
    error('No files matching C0_*.mat were found in %s', dataFolder);
end

% Preallocation
C0_mean      = nan(nCase, 1);   % Mean initial curvature from file name
dP_pert      = nan(nCase, 1);   % Initial pressure differential at perturbed cells
Vloss_pert   = nan(nCase, 1);   % Final volume change at perturbed cells

dP_eachCell    = nan(nCase, 1);
Vloss_eachCell = nan(nCase, 1);

% Loop over all curvature datasets
for kk = 1:nCase
    
    fname = fileList(kk).name;
    fpath = fullfile(fileList(kk).folder, fname);
    
    % Parse curvature value from file name
    % Example: C0_0.90667.mat --> 0.90667
    token = regexp(fname, 'C0_([0-9.]+)\.mat', 'tokens', 'once');
    if ~isempty(token)
        C0_mean(kk) = str2double(token{1});
    end
    
    % Load simulation result
    load(fpath);
    
    % Extract simulation outputs
    V_all  = Geometry{1};   % Cell volume
    P_all  = Pressure{1};   % Intracellular hydraulic pressure
    PE_all = Pressure{3};   % Interstitial hydraulic pressure
    
    % Fractional volume change relative to the initial state
    dV_all = (V_all - V_all(1,:)) ./ V_all(1,:);
    
    % Final volume change at the perturbed cell
    Vloss_cells = dV_all(end, NumPerturbed(1));
    Vloss_eachCell(kk,:) = Vloss_cells;
    Vloss_pert(kk) = mean(Vloss_cells);
    
    % Initial hydraulic pressure differential across the perturbed cell
    dP_cells = (P_all(1, NumPerturbed(1)) - PE_all(1, NumPerturbed(1))) * Pr;
    dP_eachCell(kk,:) = dP_cells;
    dP_pert(kk) = dP_cells;
end

% Sort datasets by mean initial curvature for cleaner plotting
[C0_mean_sorted, idxSort] = sort(C0_mean);
dP_sorted          = dP_pert(idxSort);
Vloss_sorted       = Vloss_pert(idxSort);
dP_each_sorted     = dP_eachCell(idxSort,:);
Vloss_each_sorted  = Vloss_eachCell(idxSort,:);

% Scatter plot: volume loss versus pressure differential
figure; hold on;

scatter(dP_sorted, 100 * Vloss_sorted, 70, C0_mean_sorted / Lr, ...
    'filled', 'MarkerFaceAlpha', 0.85, 'MarkerEdgeAlpha', 0.85);

% Optional connecting line to guide the eye
plot(dP_sorted, 100 * Vloss_sorted, '-', ...
    'Color', [0.4 0.4 0.4], 'LineWidth', 2);

xlabel('P - P^{(e)} (kPa)');
ylabel('Volume Change (%)');

set(gca, 'FontSize', 24, 'FontName', 'Arial', 'LineWidth', 1);
box on;
axis square;

cb = colorbar;
cb.Label.String = 'Mean initial curvature';
cb.FontSize = 16;
colormap(parula);