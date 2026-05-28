% -------------------------------------------------------------------------
% Figures 4, S4, and 5A–B
% Effect of initial curvature on local dynamics and response asymmetry
% -------------------------------------------------------------------------

%% Figure 4C–D: Low initial curvature (R = 2 R_t)
clc; clear;

% Load simulation data (low-curvature configuration)
load('Data_for_plot\LowCurv_R=1.mat');

% Font settings
FontSetup.FN     = 'Arial';
FontSetup.FS_txt = 22;
FontSetup.FS_cb  = 20;

% -------------------------------------------------------------------------
% Select spatial window around perturbation
% -------------------------------------------------------------------------
Idx_select_cell = NumPerturbed(1)-10 : NumPerturbed(2)+10;
Idx_select_mem  = [Idx_select_cell, Idx_select_cell(end)+1];

% Select time window (15 min before to 45 min after perturbation)
n_time_plot = 1000;
plot_time_points = linspace(timeperturb_hour - 15/tr/60, ...
                            timeperturb_hour + 45/tr/60, ...
                            n_time_plot);

temp = abs(plot_time_points - tdata');
[~, plot_time_idx] = min(temp);

% Extract local variables
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

% Plot intracellular dynamical response
Plot1_dynamical_curve_representation_Intra( ...
    tdata(plot_time_idx) - timeperturb_hour, R_all, V_all, T_all, ...
    P_all, Pi_all, PE_all, PiE_all, ...
    J_all, J_loss_all, Js_all, ...
    Pr, tr, Lr, Lseg, FontSetup);

Fig_Adjustment

% Plot initial and final configurations
figure;
plot_config_only(length(Idx_select_cell)*Lseg, length(Idx_select_cell), ...
    R_all(1,:), Rt, P_all(1,:) - PE_all(1,:), FontSetup);
clim([0,4]);

figure;
plot_config_only(length(Idx_select_cell)*Lseg, length(Idx_select_cell), ...
    R_all(end,:), Rt, P_all(end,:) - PE_all(end,:), FontSetup);
clim([0,4]);

%% Video S7: Low initial curvature ablation dynamics (R = 1.25 R_t)
% Generate videos for local zoom-in dynamics
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

FileName1 = 'LowCurv_20_cell_zoomin';
f = figure;
f.Color = 'w';
% f.Position = [600,300,800,600];

VideoWriting(length(Idx_select_cell)*Lseg, length(Idx_select_cell), Rt, ...
             tdata(plot_time_idx) - timeperturb_hour, R_select, ...
             P_select, Pi_select, PE_select, PiE_select, J_loss_select, ...
             Pr, tr, FileName1, 1, FontSetup);   % vertical colorbar

%% Figure 4G–H: High initial curvature (R = 1.25 R_t)
clc; clear;

% Load base-value dataset (higher curvature case)
load('Data_for_plot\Base_Value.mat');

FontSetup.FN     = 'Arial';
FontSetup.FS_txt = 22;
FontSetup.FS_cb  = 20;

% Same spatial and temporal selection as above
Idx_select_cell = NumPerturbed(1)-10 : NumPerturbed(end)+10;
Idx_select_mem  = [Idx_select_cell, Idx_select_cell(end)+1];

n_time_plot = 1000;
plot_time_points = linspace(timeperturb_hour - 15/tr/60, ...
                            timeperturb_hour + 45/tr/60, ...
                            n_time_plot);

temp = abs(plot_time_points - tdata');
[~, plot_time_idx] = min(temp);

% Extract local variables
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

% Plot intracellular dynamics
Plot1_dynamical_curve_representation_Intra( ...
    tdata(plot_time_idx) - timeperturb_hour, R_all, V_all, T_all, ...
    P_all, Pi_all, PE_all, PiE_all, ...
    J_all, J_loss_all, Js_all, ...
    Pr, tr, Lr, Lseg, FontSetup);

Fig_Adjustment

% Plot initial and final configurations
figure;
plot_config_only(length(Idx_select_cell)*Lseg, length(Idx_select_cell), ...
    R_all(1,:), Rt, P_all(1,:) - PE_all(1,:), FontSetup);
clim([0,22]);

figure;
plot_config_only(length(Idx_select_cell)*Lseg, length(Idx_select_cell), ...
    R_all(end,:), Rt, P_all(end,:) - PE_all(end,:), FontSetup);
clim([0,22]);

%% Figure 5A: Response asymmetry vs. mean initial curvature
clc; clear;

LOC = "Data_for_plot\Local_Varying_mean_initial_curvature\";

FontSetup.FN     = 'DejaVu Sans';
FontSetup.FS_txt = 24;
FontSetup.FS_cb  = 16;

Rt = 0.5;
C0_all = [linspace(-1,-0.0667,11), linspace(0.0667,1,11)];

for kk = 1:length(C0_all)
    
    load(strcat(LOC,'C0_',num2str(C0_all(kk)),'.mat'));
    
    % Extract curvature and compute change
    C_all  = 1 ./ Geometry{2};
    dC_all = C_all - C_all(1,:);
    
    % Absolute curvature magnitude
    Initmean_Abs  = abs(C_all);
    Changemean_Abs = abs(dC_all);
    
    % Define left/right regions relative to perturbation
    idx_right = NumPerturbed(end):NumPerturbed(end)+10;
    idx_left  = NumPerturbed(1)-10:NumPerturbed(1);
    
    % Initial curvature asymmetry index
    CAI_init(kk) = (mean(Initmean_Abs(1,idx_right)) - mean(Initmean_Abs(1,idx_left))) / ...
                   (mean(Initmean_Abs(1,idx_right)) + mean(Initmean_Abs(1,idx_left)));
    
    % Response asymmetry index
    dCAI(kk) = (mean(Changemean_Abs(end,idx_right)) - mean(Changemean_Abs(end,idx_left))) / ...
               (mean(Changemean_Abs(end,idx_right)) + mean(Changemean_Abs(end,idx_left)));
end

% Plot relationship
figure;
scatter(C0_all/Lr, dCAI, 'filled');
xlabel('Mean curvature (\mum^{-1})');
ylabel('Response Asymmetry (a.u.)');
set(gca, 'FontSize', 24, 'FontName', 'Arial', 'LineWidth', 1);
box on;
Fig_Adjustment

%% Figure 5B: Response asymmetry vs. initial curvature asymmetry
clc; clear;

C0_left_all2 = linspace(1.05, 1.95, 12);

% ---------------- Positive curvature gradient ----------------
for kk = 1:length(C0_left_all2)
    load(strcat('Data_for_plot\Local_Varying_initial_curvature_assymetry\C0_left_', ...
        num2str(C0_left_all2(kk)), '.mat'));
    
    C_all  = 1 ./ Geometry{2};
    dC_all = C_all - C_all(1,:);
    
    Initmean_Abs   = abs(C_all);
    Changemean_Abs = abs(dC_all);
    
    Left_idx  = NumPerturbed(1)-10 : NumPerturbed(1);
    Right_idx = NumPerturbed(end) : NumPerturbed(end)+10;
    
    CAI_init_Pos(kk) = (mean(Initmean_Abs(1,Right_idx)) - mean(Initmean_Abs(1,Left_idx))) / ...
                       (mean(Initmean_Abs(1,Right_idx)) + mean(Initmean_Abs(1,Left_idx)));
    
    dCAI_Pos(kk) = (mean(Changemean_Abs(end,Right_idx)) - mean(Changemean_Abs(end,Left_idx))) / ...
                   (mean(Changemean_Abs(end,Right_idx)) + mean(Changemean_Abs(end,Left_idx)));
end

% ---------------- Negative curvature gradient ----------------
for kk = 1:length(C0_left_all2)
    load(strcat('Data_for_plot\Local_Varying_initial_curvature_assymetry\Neg2_C0_left_', ...
        num2str(C0_left_all2(kk)), '.mat'));
    
    C_all  = 1 ./ Geometry{2};
    dC_all = C_all - C_all(1,:);
    
    Initmean_Abs   = abs(C_all);
    Changemean_Abs = abs(dC_all);
    
    Left_idx  = NumPerturbed(1)-10 : NumPerturbed(1);
    Right_idx = NumPerturbed(end) : NumPerturbed(end)+10;
    
    CAI_init_Neg(kk) = (mean(Initmean_Abs(1,Right_idx)) - mean(Initmean_Abs(1,Left_idx))) / ...
                       (mean(Initmean_Abs(1,Right_idx)) + mean(Initmean_Abs(1,Left_idx)));
    
    dCAI_Neg(kk) = (mean(Changemean_Abs(end,Right_idx)) - mean(Changemean_Abs(end,Left_idx))) / ...
                   (mean(Changemean_Abs(end,Right_idx)) + mean(Changemean_Abs(end,Left_idx)));
end

% ---------------- Opposing curvature configuration ----------------
for kk = 1:length(C0_left_all2)
    load(strcat('Data_for_plot\Local_Varying_initial_curvature_assymetry\Neg_C0_left_', ...
        num2str(C0_left_all2(kk)), '.mat'));
    
    C_all  = 1 ./ Geometry{2};
    dC_all = C_all - C_all(1,:);
    
    Initmean_Abs   = abs(C_all);
    Changemean_Abs = abs(dC_all);
    
    Left_idx  = NumPerturbed(1)-10 : NumPerturbed(1);
    Right_idx = NumPerturbed(end) : NumPerturbed(end)+10;
    
    CAI_init_Ops(kk) = (mean(Initmean_Abs(1,Right_idx)) - mean(Initmean_Abs(1,Left_idx))) / ...
                       (mean(Initmean_Abs(1,Right_idx)) + mean(Initmean_Abs(1,Left_idx)));
    
    dCAI_Ops(kk) = (mean(Changemean_Abs(end,Right_idx)) - mean(Changemean_Abs(end,Left_idx))) / ...
                   (mean(Changemean_Abs(end,Right_idx)) + mean(Changemean_Abs(end,Left_idx)));
end

% Plot comparison
figure;
hold on;
scatter(CAI_init_Pos, dCAI_Pos, 'filled');
scatter(CAI_init_Neg, dCAI_Neg, 'filled');
scatter(CAI_init_Ops, dCAI_Ops, 'filled');

xlabel('Initial curvature asymmetry (a.u.)');
ylabel('Response asymmetry (a.u.)');
xlim([-0.3 0.3]);

set(gca, 'FontSize', 24, 'FontName', 'Arial', 'LineWidth', 1);
box on;

Fig_Adjustment