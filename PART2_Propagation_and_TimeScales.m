% -------------------------------------------------------------------------
% Figure 3M
% Interstitial pressure kymograph and wavefront tracking
% -------------------------------------------------------------------------
clc; clear;

% These simulation results are generated from the base-value dataset
load('Data_for_plot\Base_Value.mat');
load('Data_for_plot\Data_Parsa.mat');

% Font settings
FontSetup.FN     = 'Arial';
FontSetup.FS_txt = 22;
FontSetup.FS_cb  = 16;

% Time rescaling factor
tr = 2;

%% Extract simulation variables
V_all    = Geometry{1};
R_all    = Geometry{2};
T_all    = Geometry{3};

P_all    = Pressure{1};
Pi_all   = Pressure{2};
PE_all   = Pressure{3};
PiE_all  = Pressure{4};

J_all      = Flux{1};
J_loss_all = Flux{2};
Js_all     = Flux{3};

% Quantify propagation from volume and curvature dynamics
[t_dec, d_dec] = quantify_prop(tdata, V_all, 1./R_all, NumPerturbed, timeperturb_hour, Lr, tr);

%% Interstitial pressure kymograph and wavefront tracking
% This script:
%   1. computes the interstitial pressure perturbation relative to baseline
%   2. plots a global kymograph of interstitial pressure change
%   3. plots a local kymograph centered near the perturbation site
%   4. extracts the pressure wavefront from the absolute pressure change
%   5. overlays experimental tracking data on the local kymograph
%   6. optionally estimates leftward and rightward propagation speeds

%% Select time window for visualization
n_time_plot = 100;
plot_time_points = linspace(timeperturb_hour - 3/tr/60, ...
                            timeperturb_hour + 11/tr/60, ...
                            n_time_plot);

temp = abs(plot_time_points(:) - tdata(:)');
[~, plot_time_idx] = min(temp, [], 2);
plot_time_idx = plot_time_idx(:);

%% Define spatial window around the perturbation site
N = size(PE_all, 2);
cmap = viridis(N);

Idx_select_cell = NumPerturbed(1)-20 : NumPerturbed(end)+20;
Idx_select_cell = Idx_select_cell(Idx_select_cell >= 1 & Idx_select_cell <= N);

%% Define baseline and pressure perturbation
% Interstitial pressure is measured relative to the last time point before perturbation
idx0 = find(tdata <= timeperturb_hour, 1, 'last');
if isempty(idx0)
    idx0 = 1;
end

PE0_all = PE_all(idx0,:);      % Baseline interstitial pressure
dPE_all = PE_all - PE0_all;    % Pressure perturbation relative to baseline

% Wavefront detection is based on the magnitude of the pressure change
PE_track_all = abs(dPE_all);

%% Construct plotting coordinates
x_full  = (1:N) * Lseg * Lr;
x_local = Idx_select_cell * Lseg * Lr - mean(Idx_select_cell * Lseg * Lr);
y_plot  = (tdata(plot_time_idx) - timeperturb_hour) * tr * 60;

center_x_global = mean(Idx_select_cell * Lseg * Lr);

%% Extract wavefront from the full post-perturbation trajectory
idx_start = find(tdata >= timeperturb_hour, 1, 'first');
if isempty(idx_start)
    idx_start = 1;
end

t_post = tdata(idx_start:end);
t_post = t_post(:);
t_post_min = (t_post - timeperturb_hour) * tr * 60;

PE_track_post = PE_track_all(idx_start:end, :);

% Threshold parameters for wavefront detection
threshold_fraction = 0.15;   % Detection threshold relative to global amplitude
min_response_frac  = 0.05;   % Ignore locations with very weak response

global_amp = max(PE_track_post(:));
threshold_value = threshold_fraction * global_amp;

% WaveFront_PE stores [position, crossing time]
WaveFront_PE = nan(N, 2);

for j = 1:N
    trace_j = PE_track_post(:, j);
    trace_j = trace_j(:);

    % Ignore locations with negligible response
    if max(trace_j) < min_response_frac * global_amp
        continue;
    end

    % Find first threshold crossing
    idx_cross = find(trace_j >= threshold_value, 1, 'first');
    if isempty(idx_cross)
        continue;
    end

    % Linearly interpolate the threshold-crossing time
    if idx_cross == 1
        t_cross = t_post_min(1);
    else
        t1 = t_post_min(idx_cross-1);
        t2 = t_post_min(idx_cross);
        y1 = trace_j(idx_cross-1);
        y2 = trace_j(idx_cross);

        if abs(y2 - y1) < eps
            t_cross = t2;
        else
            t_cross = t1 + (threshold_value - y1) * (t2 - t1) / (y2 - y1);
        end
    end

    WaveFront_PE(j,:) = [x_full(j), t_cross];
end

% Keep only valid wavefront points
valid_wave = ~isnan(WaveFront_PE(:,2));
WaveFront_PE = WaveFront_PE(valid_wave,:);

%% Convert wavefront into local coordinates
WaveFront_PE_local = WaveFront_PE;
WaveFront_PE_local(:,1) = WaveFront_PE(:,1) - mean(Idx_select_cell * Lseg * Lr);

%% Split left and right branches
left_idx_full  = WaveFront_PE(:,1) < center_x_global;
right_idx_full = WaveFront_PE(:,1) > center_x_global;

WaveFront_PE_left_full  = sortrows(WaveFront_PE(left_idx_full,:), 1);
WaveFront_PE_right_full = sortrows(WaveFront_PE(right_idx_full,:), 1);

left_idx_local  = WaveFront_PE_local(:,1) < 0;
right_idx_local = WaveFront_PE_local(:,1) > 0;

WaveFront_PE_left_local  = sortrows(WaveFront_PE_local(left_idx_local,:), 1);
WaveFront_PE_right_local = sortrows(WaveFront_PE_local(right_idx_local,:), 1);

%% Build wavefront curves for plotting: global coordinates
WaveFront_PE_left_full_plot  = [];
WaveFront_PE_right_full_plot = [];
WaveFront_PE_full_plot       = [];

if ~isempty(WaveFront_PE_left_full)
    valid_left = WaveFront_PE_left_full(:,2) >= min(y_plot) & ...
                 WaveFront_PE_left_full(:,2) <= max(y_plot);
    WaveFront_PE_left_full_plot = WaveFront_PE_left_full(valid_left,:);
end

if ~isempty(WaveFront_PE_right_full)
    valid_right = WaveFront_PE_right_full(:,2) >= min(y_plot) & ...
                  WaveFront_PE_right_full(:,2) <= max(y_plot);
    WaveFront_PE_right_full_plot = WaveFront_PE_right_full(valid_right,:);
end

if ~isempty(WaveFront_PE_left_full_plot) && ~isempty(WaveFront_PE_right_full_plot)
    WaveFront_PE_full_plot = [WaveFront_PE_left_full_plot; WaveFront_PE_right_full_plot];
elseif ~isempty(WaveFront_PE_left_full_plot)
    WaveFront_PE_full_plot = WaveFront_PE_left_full_plot;
elseif ~isempty(WaveFront_PE_right_full_plot)
    WaveFront_PE_full_plot = WaveFront_PE_right_full_plot;
end

%% Build wavefront curves for plotting: local coordinates
WaveFront_PE_left_local_plot  = [];
WaveFront_PE_right_local_plot = [];
WaveFront_PE_local_plot       = [];

if ~isempty(WaveFront_PE_left_local)
    valid_left = WaveFront_PE_left_local(:,1) >= min(x_local) & ...
                 WaveFront_PE_left_local(:,1) <= max(x_local) & ...
                 WaveFront_PE_left_local(:,2) >= min(y_plot) & ...
                 WaveFront_PE_left_local(:,2) <= max(y_plot);
    WaveFront_PE_left_local_plot = WaveFront_PE_left_local(valid_left,:);
end

if ~isempty(WaveFront_PE_right_local)
    valid_right = WaveFront_PE_right_local(:,1) >= min(x_local) & ...
                  WaveFront_PE_right_local(:,1) <= max(x_local) & ...
                  WaveFront_PE_right_local(:,2) >= min(y_plot) & ...
                  WaveFront_PE_right_local(:,2) <= max(y_plot);
    WaveFront_PE_right_local_plot = WaveFront_PE_right_local(valid_right,:);
end

if ~isempty(WaveFront_PE_left_local_plot) && ~isempty(WaveFront_PE_right_local_plot)
    WaveFront_PE_local_plot = [WaveFront_PE_left_local_plot; WaveFront_PE_right_local_plot];
elseif ~isempty(WaveFront_PE_left_local_plot)
    WaveFront_PE_local_plot = WaveFront_PE_left_local_plot;
elseif ~isempty(WaveFront_PE_right_local_plot)
    WaveFront_PE_local_plot = WaveFront_PE_right_local_plot;
end

%% Global kymograph of interstitial pressure perturbation
figure;
PE_global_plot = dPE_all(plot_time_idx, :);

% Downsample for cleaner vector-style rendering
n_time_ds  = 200;
n_space_ds = 150;

idx_t = round(linspace(1, size(PE_global_plot,1), n_time_ds));
idx_t = unique(idx_t);

idx_x = round(linspace(1, size(PE_global_plot,2), n_space_ds));
idx_x = unique(idx_x);

PE_ds = PE_global_plot(idx_t, idx_x);
x_ds  = x_full(idx_x);
y_ds  = y_plot(idx_t);

h = pcolor(x_ds, y_ds, PE_ds);
set(h, 'EdgeColor', 'none');
set(gca, 'YDir', 'normal');
colormap(cmap);
clim([-1 3]);

cb = colorbar;
cb.Label.String = '\Delta P_E';

xlabel('AP Position (\mum)');
ylabel('Time (min)');
title('Interstitial Pressure Perturbation \Delta P_E', 'FontWeight', 'Normal');
set(gca, 'FontSize', FontSetup.FS_txt, 'FontName', FontSetup.FN, 'LineWidth', 1);

hold on;

h_wave_full = [];
if ~isempty(WaveFront_PE_full_plot)
    h_wave_full = plot(WaveFront_PE_full_plot(:,1), WaveFront_PE_full_plot(:,2), ...
        '-', 'Color', [1 1 1], 'LineWidth', 2.5, ...
        'DisplayName', 'Interstitial pressure wavefront');
end

Fig_Adjustment

%% Local kymograph of interstitial pressure perturbation
figure;
PE_local_plot = dPE_all(plot_time_idx, Idx_select_cell);

h = pcolor(x_local, y_plot, PE_local_plot);
set(h, 'EdgeColor', 'none');
colormap(cmap);

maxAbsPE_local = max(abs(PE_local_plot(:)));
if maxAbsPE_local > 0
    clim([-maxAbsPE_local, maxAbsPE_local]);
end

cb = colorbar;
cb.Label.String = '\Delta P_E';
clim([-1 3]);

xlabel('Position (\mum)');
ylabel('Time (min)');
title('Local Interstitial Pressure Perturbation \Delta P_E', 'FontWeight', 'Normal');
set(gca, 'FontSize', FontSetup.FS_txt, 'FontName', FontSetup.FN, 'LineWidth', 1);

hold on;

h_wave = [];
% If desired, the simulated local wavefront can also be overlaid here
% if ~isempty(WaveFront_PE_local_plot)
%     h_wave = plot(WaveFront_PE_local_plot(:,1), WaveFront_PE_local_plot(:,2), ...
%         '-', 'Color', [1 1 1], 'LineWidth', 2.5, ...
%         'DisplayName', 'Interstitial pressure wavefront');
% end

h_exp = [];
if exist('HalfMaxVolTrack', 'var') && ~isempty(HalfMaxVolTrack)
    h_exp = scatter(HalfMaxVolTrack(:,1), HalfMaxVolTrack(:,2), ...
        50, 'filled', ...
        'MarkerFaceColor', [1 0.5 0], ...
        'LineWidth', 0.8, ...
        'DisplayName', 'Experiment');
end

% Add legend only for available plotted objects
if ~isempty(h_wave) && ~isempty(h_exp)
    legend([h_wave, h_exp], ...
        {'Interstitial pressure wavefront', 'Experiment'}, ...
        'Location', 'best', 'Box', 'off');
elseif ~isempty(h_wave)
    legend(h_wave, {'Interstitial pressure wavefront'}, ...
        'Location', 'best', 'Box', 'off');
elseif ~isempty(h_exp)
    legend(h_exp, {'Experiment'}, ...
        'Location', 'best', 'Box', 'off');
end

Fig_Adjustment

%% Optional: estimate leftward and rightward propagation speeds
if ~isempty(WaveFront_PE_left_local_plot)
    x_left = abs(WaveFront_PE_left_local_plot(:,1));
    t_left = WaveFront_PE_left_local_plot(:,2);

    if numel(x_left) >= 2
        p_left = polyfit(x_left, t_left, 1);   % t = a*x + b
        speed_left = 1 / p_left(1);            % um/min
        fprintf('Left PE wavefront speed = %.4f um/min\n', speed_left);
    end
end

if ~isempty(WaveFront_PE_right_local_plot)
    x_right = WaveFront_PE_right_local_plot(:,1);
    t_right = WaveFront_PE_right_local_plot(:,2);

    if numel(x_right) >= 2
        p_right = polyfit(x_right, t_right, 1);   % t = a*x + b
        speed_right = 1 / p_right(1);             % um/min
        fprintf('Right PE wavefront speed = %.4f um/min\n', speed_right);
    end
end

%The updated propagation velocity (calculate the mean velocity in the first 2 min after ablation)
v_new = (x_left(4)/t_left(4)+x_right(4)/t_right(4))/2;