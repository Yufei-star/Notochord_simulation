% -------------------------------------------------------------------------
% Preprocessing and smoothing of experimental volume data
% -------------------------------------------------------------------------

% List of variable names
load('Volume_all_dpf_exp.mat');
varNames = {'V10dpf','V2dpf','V3dpf','V4dpf','V5dpf','V6dpf','V7dpf','V8dpf','V9dpf'};

scale = 0.7279847; % um per pixel

% Smoothing parameters
method = 'movmean';   % 'sgolay' or 'movmean'
window = 40;         % must be odd for sgolay
polyOrder = 3;       % for sgolay

for i = 1:length(varNames)
    
    % Get variable
    data = eval(varNames{i});
    
    % Remove NaNs
    data = data(~any(isnan(data),2), :);
    
    % Extract
    x = data(:,1);
    V = data(:,2);
    
    % Convert x to microns
    x = x * scale;
    
    % Shift to start at 0
    x = x - min(x);
    
    % Sort
    [x_sorted, idx] = sort(x);
    V_sorted = V(idx);
    
    % ---------- SMOOTHING ----------
    switch method
        case 'sgolay'
            V_smooth = sgolayfilt(V_sorted, polyOrder, window);
        case 'movmean'
            V_smooth = smoothdata(V_sorted, 'movmean', window);
    end
    
    % Store processed data
    data_proc = [x_sorted, V_smooth];
    
    % Save
    newName = [varNames{i}, '_proc'];
    assignin('base', newName, data_proc);
    
end