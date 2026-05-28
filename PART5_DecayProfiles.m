% -------------------------------------------------------------------------
% Figure 5D and Figure S5
% Decay length scale and parameter dependence of response dynamics
% -------------------------------------------------------------------------

%% Figure 5D: Decay length scale of curvature response
clc; clear;

LOC = "Data_for_plot\Global_Varying_mean_initial_curvature\";

FontSetup.FN     = 'Arial';
FontSetup.FS_txt = 24;
FontSetup.FS_cb  = 16;

Rt = 0.5;

% Define initial curvature range through initial radius
R0_all_temp = linspace(1,2,15) * Rt;
R0_all2 = [R0_all_temp(4:end), 2:15];

for kkk = 1:length(R0_all2)
    
    load(strcat(LOC, 'R0_', num2str(R0_all2(kkk)), '.mat'));
    
    % Extract variables
    V_all = Geometry{1};
    R_all = Geometry{2};
    C_all = 1 ./ R_all;
    
    % Identify perturbation time index
    [~, idx] = min(abs(tdata - timeperturb_hour));
    
    % Compute normalized volume and curvature responses
    dVm_all = max(abs(V_all(idx:end,:) - V_all(1,:)) ./ V_all(1,:));
    dCm_all = max(abs(C_all(idx:end,:) - C_all(1,:))) / Lr;
    
    % Normalize curvature response and remove baseline offset
    dCm_all = dCm_all - dCm_all(end);
    
    % Use right-hand side for decay fitting
    dC_mag_right = dCm_all(NumPerturbed(end)+1:end);
    X = (1:length(dC_mag_right)) * Lseg * Lr;
    
    % Remove invalid or near-zero values
    valid = isfinite(X) & isfinite(dC_mag_right) & (dC_mag_right > 1e-8);
    X_fit = X(valid);
    Y_fit = dC_mag_right(valid);
    
    % Sort for fitting
    [X_fit, idx_sort] = sort(X_fit);
    Y_fit = Y_fit(idx_sort);
    
    % ---------------- Exponential decay fitting ----------------
    ft = fittype('A*exp(lambda*x)', ...
        'independent','x', ...
        'coefficients',{'A','lambda'});
    
    % Initial guesses
    A0 = Y_fit(1);
    lambda0 = -1/(max(X_fit) + eps);
    
    % Fit exponential decay
    [fitobj, ~] = fit(X_fit(:), Y_fit(:), ft, ...
        'StartPoint', [A0, lambda0], ...
        'Lower', [0, -Inf], ...
        'Upper', [Inf, 0]);
    
    % Extract decay length scale
    lambda = fitobj.lambda;
    LengthScale_all(kkk) = -1 / lambda;
end

% Plot decay length scale vs curvature
figure;
plot(1 ./ (R0_all2(2:end) * Lr), LengthScale_all(2:end), 'LineWidth', 2);

xlabel('Curvature (μm^{-1})');
ylabel('Length scale (μm)');
set(gca, 'FontSize', FontSetup.FS_txt, 'FontName', FontSetup.FN, 'LineWidth', 1);

%% Figure S5A–B: Influence of membrane permeability
clc; clear;

device_name = 'Data_for_plot\';
FontSetup.FN     = 'Arial';
FontSetup.FS_txt = 22;
FontSetup.FS_cb  = 20;

alpha_isf_scan = [0 0.001 0.01 0.1];
device = [device_name, 'alpha_isf_'];
FigNum = [1, 2];

for kk = 1:length(alpha_isf_scan)
    
    FileName = [device, num2str(alpha_isf_scan(kk)), '.mat'];
    load(FileName);
    
    % Plot spatial decay profiles
    Plot_Decay_Profile(tdata, Geometry, Pressure, Flux, ...
        timeperturb_hour, NumPerturbed, Lr, Lseg, FigNum, FontSetup);
    
    % Compute response amplitudes
    V_all = Geometry{1};
    R_all = Geometry{2};
    C_all = 1 ./ R_all;
    
    dV_all = (V_all - V_all(1,:)) ./ V_all(1,:);
    dC_all = (C_all - C_all(1,:)) / Lr;
    
    dVmax_permeability(kk) = min(dV_all(end,:));
    dCmax_permeability(kk) = max(abs(dC_all(:)));
end

% Add legend
figure(FigNum(1));
l = legend('$\tilde{\alpha} = 0$', '$\tilde{\alpha} = 10^{-3}$', ...
           '$\tilde{\alpha} = 10^{-2}$', '$\tilde{\alpha} = 10^{-1}$', ...
           'Interpreter', 'latex');
set(l, 'box', 'off', 'FontSize', 18);

Fig_Adjustment

%% Figure S5C-D: Bar plots of response amplitude vs permeability
co = get(groot, 'defaultAxesColorOrder');
c1 = co(1,:);
c2 = co(2,:);

% Volume response
figure;
b = bar(dVmax_permeability(:) * 100);
b.FaceColor = c1;
b.EdgeColor = c1;

set(gca, 'FontSize', 18, 'FontName', 'Arial');
xticks(1:length(alpha_isf_scan));
xticklabels({'0','10^{-3}','10^{-2}','10^{-1}'});
xlabel('$\tilde{\alpha}$', 'Interpreter', 'latex');
ylabel('\Delta V (%)');
axis square; box on;

% Curvature response
figure;
b = bar(dCmax_permeability(:));
b.FaceColor = c2;
b.EdgeColor = c2;

set(gca, 'FontSize', 18, 'FontName', 'Arial');
xticks(1:length(alpha_isf_scan));
xticklabels({'0','10^{-3}','10^{-2}','10^{-1}'});
xlabel('$\tilde{\alpha}$', 'Interpreter', 'latex');
ylabel('\Delta C_{max} (\mu m^{-1})');
axis square; box on;

%% Figure S5E-F: Influence of initial configuration
FNs = ["C0_0.0667.mat","C0_1.mat","Neg_C0_left_1.5409.mat"];
device = 'Data_for_plot\Profile_LengthScale\';
FigNum = [7, 8];

for kk = 1:length(FNs)
    
    load(strcat(device, FNs(kk)));
    
    % Plot decay profiles
    Plot_Decay_Profile(tdata, Geometry, Pressure, Flux, ...
        timeperturb_hour, NumPerturbed, Lr, Lseg, FigNum, FontSetup);
    
    % Compute response amplitudes
    V_all = Geometry{1};
    R_all = Geometry{2};
    C_all = 1 ./ R_all;
    
    dV_all = (V_all - V_all(1,:)) ./ V_all(1,:);
    dC_all = (C_all - C_all(1,:)) / Lr;
    
    dVmax_IC(kk) = min(dV_all(end,:));
    dCmax_IC(kk) = max(abs(dC_all(:)));
end

figure(FigNum(1));
l = legend('Low \langle C \rangle', ...
           'High \langle C \rangle', ...
           'Opposing curvature');
set(l, 'box', 'off', 'FontSize', 18);

Fig_Adjustment

% Plot representative initial configurations
for kk = 1:length(FNs)
    
    load(strcat(device, FNs(kk)));
    figure(90 + kk);
    
    % Select local window
    Idx_select_cell = NumPerturbed(1)-10 : NumPerturbed(end)+10;
    Idx_select_mem  = [Idx_select_cell, Idx_select_cell(end)+1];
    
    % Time selection
    n_time_plot = 1000;
    plot_time_points = linspace(timeperturb_hour - 15/tr/60, ...
                                timeperturb_hour + 45/tr/60, ...
                                n_time_plot);
    
    temp = abs(plot_time_points - tdata');
    [~, plot_time_idx] = min(temp);
    
    % Extract local geometry
    R_all = Geometry{2}(plot_time_idx, Idx_select_mem);
    P_all = Pressure{1}(plot_time_idx, Idx_select_cell);
    PE_all   = Pressure{3}(plot_time_idx, Idx_select_cell);
    
    % Plot configuration
    plot_config_only(length(Idx_select_cell)*Lseg, length(Idx_select_cell), ...
        R_all(1,:), Rt, P_all(1,:)-PE_all(1,:), FontSetup);
end

%% Figure S5G–H: Bar plots for different initial configurations
co = get(groot, 'defaultAxesColorOrder');
c1 = co(1,:);
c2 = co(2,:);

labels_IC = {'Low \langle C \rangle', 'High \langle C \rangle', 'Opposing curvature'};

% Volume response
figure;
b = bar(dVmax_IC(:) * 100);
b.FaceColor = c1;
b.EdgeColor = c1;

set(gca, 'FontSize', 18, 'FontName', 'Arial');
xticks(1:3);
xticklabels(labels_IC);
ylabel('\Delta V (%)');
axis square; box on;

% Curvature response
figure;
b = bar(dCmax_IC(:));
b.FaceColor = c2;
b.EdgeColor = c2;

set(gca, 'FontSize', 18, 'FontName', 'Arial');
xticks(1:3);
xticklabels(labels_IC);
ylabel('\Delta C_{max} (\mu m^{-1})');
axis square; box on;