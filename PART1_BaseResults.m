% -------------------------------------------------------------------------
% Figures 3, S3, and Videos S5-S6
% Base-case laser ablation simulation: local dynamics, kymographs, and videos
% -------------------------------------------------------------------------
clc; clear;
load('Data_for_plot\Base_Value.mat');

% Font settings for all figures
FontSetup.FN     = 'Arial';
FontSetup.FS_txt = 22;
FontSetup.FS_cb  = 20;

%% Figure 3H,I and Figure S3F
% Base-case local dynamics near the puncture site:
% volume, curvature, pressure, flux, and membrane tension

% Select a 20-cell window centered around the perturbed region
Idx_select_cell = NumPerturbed(1)-10 : NumPerturbed(end)+10;
Idx_select_mem  = [Idx_select_cell, Idx_select_cell(end)+1];

% Select time points spanning 15 min before to 45 min after perturbation
n_time_plot = 1000;
plot_time_points = linspace(timeperturb_hour - 15/tr/60, ...
                            timeperturb_hour + 45/tr/60, ...
                            n_time_plot);
temp = abs(plot_time_points - tdata');
[~, plot_time_idx] = min(temp);

% Extract local variables within the selected spatial window
V_all    = Geometry{1}(plot_time_idx, Idx_select_cell);
R_all    = Geometry{2}(plot_time_idx, Idx_select_mem);
T_all    = Geometry{3}(plot_time_idx, Idx_select_mem);

P_all    = Pressure{1}(plot_time_idx, Idx_select_cell);
Pi_all   = Pressure{2}(plot_time_idx, Idx_select_cell);
PE_all   = Pressure{3}(plot_time_idx, Idx_select_cell);
PiE_all  = Pressure{4}(plot_time_idx, Idx_select_cell);

J_all      = Flux{1}(plot_time_idx, Idx_select_cell);
J_loss_all = Flux{2}(plot_time_idx, Idx_select_cell);
Js_all     = Flux{3}(plot_time_idx, Idx_select_cell);

% Plot intracellular and local dynamical responses as line plots
Plot1_dynamical_curve_representation_Intra( ...
    tdata(plot_time_idx) - timeperturb_hour, R_all, V_all, T_all, ...
    P_all, Pi_all, PE_all, PiE_all, ...
    J_all, J_loss_all, Js_all, ...
    Pr, tr, Lr, Lseg, FontSetup);

% Unify x-axis limits across all generated axes
figs = findall(0, 'Type', 'figure');
for f = figs'
    ax = findall(f, 'Type', 'axes', '-not', 'Tag', 'Colorbar', '-not', 'Tag', 'legend');
    for a = ax'
        xlim(a, [-2 10]);
    end
end

Fig_Adjustment

%% Figure 3K,L
% Local intracellular and interstitial dynamics shown as 2D kymographs

% Select the same 20-cell window around the perturbed region
Idx_select_cell = NumPerturbed(1)-10 : NumPerturbed(2)+10;
Idx_select_mem  = [Idx_select_cell, Idx_select_cell(end)+1];

% Select time points spanning 3 min before to 10 min after perturbation
n_time_plot = 100;
plot_time_points = linspace(timeperturb_hour - 3/tr/60, ...
                            timeperturb_hour + 10/tr/60, ...
                            n_time_plot);
temp = abs(plot_time_points - tdata');
[~, plot_time_idx] = min(temp);

% Extract local variables for kymograph plotting
V_select    = Geometry{1}(plot_time_idx, Idx_select_cell);
R_select    = Geometry{2}(plot_time_idx, Idx_select_mem);
T_select    = Geometry{3}(plot_time_idx, Idx_select_mem);

P_select    = Pressure{1}(plot_time_idx, Idx_select_cell);
Pi_select   = Pressure{2}(plot_time_idx, Idx_select_cell);
PE_select   = Pressure{3}(plot_time_idx, Idx_select_cell);
PiE_select  = Pressure{4}(plot_time_idx, Idx_select_cell);

J_select      = Flux{1}(plot_time_idx, Idx_select_cell);
J_loss_select = Flux{2}(plot_time_idx, Idx_select_cell);
Js_select     = Flux{3}(plot_time_idx, Idx_select_cell);

% Also extract full-tissue variables for whole-system visualization
V_all    = Geometry{1}(plot_time_idx, :);
R_all    = Geometry{2}(plot_time_idx, :);
T_all    = Geometry{3}(plot_time_idx, :);

P_all    = Pressure{1}(plot_time_idx, :);
Pi_all   = Pressure{2}(plot_time_idx, :);
PE_all   = Pressure{3}(plot_time_idx, :);
PiE_all  = Pressure{4}(plot_time_idx, :);

J_all      = Flux{1}(plot_time_idx, :);
J_loss_all = Flux{2}(plot_time_idx, :);
Js_all     = Flux{3}(plot_time_idx, :);

% Plot local intra- and interstitial dynamics as line plots and kymographs
Plot2_2D_kymograph_Intra_and_Interstitial( ...
    tdata(plot_time_idx) - timeperturb_hour, R_select, V_select, ...
    P_select, Pi_select, PE_select, ...
    Pr, tr, Lr, Lseg, FontSetup);

Fig_Adjustment

%% Videos S5-S6
% Generate videos for local zoom-in and whole-tissue dynamics

% Video S5: 20-cell zoom-in around the perturbation site
FileName1 = 'BaseValue_20_cell_zoomin_Arial';
f = figure;
f.Color = 'w';
% f.Position = [600,300,800,600];

VideoWriting(length(Idx_select_cell)*Lseg, length(Idx_select_cell), Rt, ...
             tdata(plot_time_idx) - timeperturb_hour, R_select, ...
             P_select, Pi_select, PE_select, PiE_select, J_loss_select, ...
             Pr, tr, FileName1, 1, FontSetup);   % vertical colorbar

% Video S6: whole-tissue dynamics
FileName2 = 'BaseValue_all_cells_Arial';
f = figure;
f.Color = 'w';
f.Position = [600,300,800,600];

VideoWriting(N*Lseg, N, Rt, ...
             tdata(plot_time_idx) - timeperturb_hour, R_all, ...
             P_all, Pi_all, PE_all, PiE_all, J_loss_all, ...
             Pr, tr, FileName2, 0, FontSetup);   % horizontal colorbar