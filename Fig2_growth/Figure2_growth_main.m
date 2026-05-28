% -------------------------------------------------------------------------
% Run the model and generate simulation data
% -------------------------------------------------------------------------

clc; clear;

% Plot settings
FN = 'Arial';
FontSize_text = 22;

% Base parameter values
N = 20;                         % Number of cells
n_mesh = 1:N;
Lseg = 0.6*ones(1,N) + 4e-4*(n_mesh-1).*(N-n_mesh+3);   % Initial cell lengths

k = 0.5;                        % Membrane elastic modulus
alpha_isf = 1e-2 * ones(1, N);  % Permeability of cell membrane (lateral surface, interface between intersitial fluid)
alpha_cell = zeros(1, N);       % Permeability between neighboring cells (left membrane)
alphaw_ext = 2;                 % Water transport in interstitial fluid (Darcy-type transport)

gamma_perturbed = 4.5e-3 + 2e-3*(linspace(0,1,N)).^0.5; % Spatially varying ion import rate
vL = 1.1e-2*(1 + 0.34*(linspace(0,1,N)).^(0.5));        % Spatially varying cell elongation rate

Pr_in = 10;                     % intracellular pressure at the tail cell (i=20)
time = 12;                      % Total simulation time (dimensionless)
timeperturb_hour = 0;           % Perturbation time point (h)

% Initial geometric configuration
Rt = 0.5 * 1.25;
R0 = 2.5 * Rt;             % Initial radius
R = R0 * ones(1, N+1);          % Boundary radii
R(1) = -R0;                     % Left boundary
kp = 1 ./ R;

% Characteristic scales
Lr = 2 * 20;    % Length scale = 2Rt (um)
Pr = 1;         % Pressure scale (kPa)
tr = 2;     % Time scale (h)

% Run simulation
CaseName = 'Growth pattern';
[tdata, Geometry, Pressure, Flux] = Notochord_main_pkg_long( ...
    N, Lseg, alpha_isf, gamma_perturbed, vL, ...
    alphaw_ext, alpha_cell, k, Rt, R, Pr_in, ...
    time, timeperturb_hour, CaseName);

% Select time points for plotting
n_time_plot = 100;
plot_time_points = linspace(0, time, n_time_plot);
temp = abs(plot_time_points - tdata');
[~, plot_time_idx] = min(temp);

% Extract simulation outputs
V_all    = Geometry{1}(plot_time_idx,:);
R_all    = Geometry{2}(plot_time_idx,:);
T_all    = Geometry{3}(plot_time_idx,:);
Lseg_all = Geometry{4}(plot_time_idx,:);

P_all   = Pressure{1}(plot_time_idx,:);
Pi_all  = Pressure{2}(plot_time_idx,:);
PE_all  = Pressure{3}(plot_time_idx,:);
PiE_all = Pressure{4}(plot_time_idx,:);

J_all      = Flux{1}(plot_time_idx,:);
J_loss_all = Flux{2}(plot_time_idx,:);
Js_all     = Flux{3}(plot_time_idx,:);

C_all = 1 ./ R_all;

% Compute derived quantities for growth analysis
GR = log(V_all(end,:) ./ V_all(1,:));        % Volume growth rate
dC = C_all(end,:) - C_all(1,:);              % Curvature change

xmesh = cumsum(Lseg) * Lr;
XLIM = [0 - Lseg_all(1,1), sum(Lseg_all(1,:)) * 1.275];

%% Plot growth rate, volume profile, and curvature profile from 2 dpf to 3 dpf

% -------------------------------------------------------------------------
% Fig. 2F: Volume profiles from simulation and experiment at 48 hpf and 72 hpf
% -------------------------------------------------------------------------
figure; hold on;

% Load experimental volume data
load('Volume_all_dpf_exp_proc_smooth_movmean40.mat', 'V2dpf_proc', 'V3dpf_proc');

% Simulation / theory profiles
h_sim1 = plot((cumsum(Lseg_all(1,2:end-1)) - Lseg_all(1,2)) * Lr + 600, ...
              1e-4 * V_all(1,2:end-1) * Lr^3, ...
              'linewidth', 2);

h_sim3 = plot((cumsum(Lseg_all(end,2:end-1)) - Lseg_all(end,2)) * Lr + 600, ...
              1e-4 * V_all(end,2:end-1) * Lr^3, ...
              'linewidth', 2);

% Colors from simulation lines
c48 = h_sim1.Color;
c72 = h_sim3.Color;

% Experimental data: 2 dpf = 48 hpf
xexp2 = V2dpf_proc(:,1);
Vexp2 = V2dpf_proc(:,2);

[xexp2, ia2] = unique(xexp2, 'stable');
Vexp2 = Vexp2(ia2);

idx2 = (xexp2 >= 600) & (xexp2 <= 1000);

s2 = scatter(xexp2(idx2) - min(xexp2(idx2)) + 600, ...
             1e-4 * Vexp2(idx2), ...
             40, 'o', ...
             'MarkerFaceColor', c48, ...
             'MarkerEdgeColor', c48, ...
             'MarkerFaceAlpha', 0.35, ...
             'MarkerEdgeAlpha', 0.35, ...
             'LineWidth', 1.0);

% Experimental data: 3 dpf = 72 hpf
xexp3 = V3dpf_proc(:,1);
Vexp3 = V3dpf_proc(:,2);

[xexp3, ia3] = unique(xexp3, 'stable');
Vexp3 = Vexp3(ia3);

idx3 = (xexp3 >= 700) & (xexp3 <= 1250);

s3 = scatter(xexp3(idx3) - min(xexp3(idx3)) + 600, ...
             1e-4 * Vexp3(idx3), ...
             40, 'o', ...
             'MarkerFaceColor', c72, ...
             'MarkerEdgeColor', c72, ...
             'MarkerFaceAlpha', 0.35, ...
             'MarkerEdgeAlpha', 0.35, ...
             'LineWidth', 1.0);

xlabel('AP Position (\mum)');
ylabel('Cell Volume (\times 10^4 \mum^3)');

l1 = legend([h_sim1, s2, h_sim3, s3], ...
            {'48 hpf (theory)', '48 hpf (exp.)', ...
             '72 hpf (theory)', '72 hpf (exp.)'}, ...
            'Location', 'southeast', ...
            'NumColumns', 1);
set(l1, 'box', 'off');

set(gca, 'FontSize', FontSize_text, 'Fontname', FN, 'linewidth', 1);

% -------------------------------------------------------------------------
% Fig. 2G: Curvature profiles at the beginning and end of development
% -------------------------------------------------------------------------
figure;
plot(xmesh(1:end-1) + 600, 1e2 * C_all(1,2:end-1) / Lr, 'linewidth', 2); hold on;
plot(xmesh(1:end-1) + 600, 1e2 * C_all(end,2:end-1) / Lr, 'linewidth', 2);
yline(0, '--', 'LineWidth', 1.5);   % Zero-curvature reference

xlabel('AP Position (μm)');
ylabel('Curvature (10^{-2}\times μm^{-1})');
xlim([min(xmesh(1:end-1)), max(xmesh(1:end-1))]);
ylim(1e2 * [-0.04 0.04]);

l2 = legend('24 hpf', '48 hpf', 'Location', 'southeast', 'NumColumns', 2);
set(l2, 'box', 'off');

set(gca, 'FontSize', FontSize_text, 'Fontname', FN, 'linewidth', 1);

% -------------------------------------------------------------------------
% Fig. 2H: Volume growth rate and curvature change along the AP axis
% -------------------------------------------------------------------------
figure;
yyaxis left
plot(xmesh(1:end-1) + 600, GR(1:end-1), 'linewidth', 2); hold on;
xlabel('AP Position (μm)');
ylabel('Growth Rate');

yyaxis right
plot(xmesh(1:end-1) + 600, 1e2 * dC(2:end-1) / Lr, 'linewidth', 2);
ylabel('\DeltaC (\times10^{-2}μm^{-1})');
xlim([min(xmesh(1:end-1)) max(xmesh(1:end-1))]);

set(gca, 'FontSize', FontSize_text, 'Fontname', FN, 'linewidth', 1);

% Fig. 2D-E: Hydraulic and osmotic pressure kymographs
cmap = viridis(N);

% -------------------------------------------------------------------------
% Fig. 2D: Hydraulic pressure difference kymograph
% -------------------------------------------------------------------------
figure;
h2 = pcolor(cumsum(Lseg) * Lr + 600, tdata(plot_time_idx) * tr, (P_all - PE_all) * Pr);
shading interp;
colormap(cmap);
colorbar;
set(h2, 'Edgecolor', 'None');

xlabel('AP Position (μm)');
ylabel('Time (h)');
title('P - P_0 (kPa)', 'FontWeight', 'Normal');

set(gca, 'FontSize', FontSize_text, 'Fontname', FN, 'linewidth', 1);

% -------------------------------------------------------------------------
% Fig. 2E: Osmotic pressure difference kymograph
% -------------------------------------------------------------------------
figure;
h4 = pcolor(cumsum(Lseg(2:N)) * Lr + 600, tdata(plot_time_idx) * tr, (Pi_all(:,2:N) - PiE_all(:,2:N)) * Pr);
shading interp;
colormap(cmap);
colorbar;
set(h4, 'Edgecolor', 'None');

xlabel('AP Position (μm)');
ylabel('Time (h)');
title('\Pi - \Pi_0 (kPa)', 'FontWeight', 'Normal');

set(gca, 'FontSize', FontSize_text, 'Fontname', FN, 'linewidth', 1);

% Figure adjustment for all figures above
Fig_Adjustment

% -------------------------------------------------------------------------
% Fig. 2C: Configuration plots for
% -------------------------------------------------------------------------
FontSetup.FN = 'Arial';
FontSetup.FS_cb = 18;
FontSetup.FS_txt = 16;

figure;
plot_config_only(Lseg_all(1,:), N, R_all(1,:), Rt, P_all(1,:) - PE_all(1,:), XLIM, FontSetup);
clim([10.8 11.65] - 1);

figure;
plot_config_only(Lseg_all(end,:), N, R_all(end,:), Rt, P_all(end,:) - PE_all(end,:), XLIM, FontSetup);
clim([10.8 11.65] - 1);

%% Video S2: Growth simulation with hydraulic and osmotic pressure dynamics
FontSetup.FN = 'Arial';
FontSetup.FS_cb = 12;
FontSetup.FS_txt = 16;

FileName1 = 'Growth pattern';

f = figure;
f.Color = 'w';
f.Position = [600,300,600,400];

VideoWriting_Growth(Lseg, N, Rt, tdata(plot_time_idx) - timeperturb_hour, R_all, Lseg_all, ...
                    P_all - PE_all, Pi_all - PiE_all, PE_all, PiE_all, J_loss_all, ...
                    Pr, tr, FileName1, 1, FontSetup)