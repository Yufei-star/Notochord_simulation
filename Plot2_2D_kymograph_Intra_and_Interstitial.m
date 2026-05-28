function [] = Plot2_2D_kymograph_Intra_and_Interstitial( ...
    tdata,R_all,V_all,P_all,Pi_all,PE_all,Pr,tr,Lr,Lseg,FontSetup)

% =========================
% Settings
% =========================
n_time_plot = 300;   % reduce from 1000 to 300 (tune: 200~400 is usually enough)

video_time_points = linspace(min(tdata), max(tdata), n_time_plot);
temp = abs(video_time_points - tdata');
[~, plot_time_idx] = min(temp);

[~, N] = size(V_all);
cmap = viridis(256);   % use fixed colormap resolution instead of viridis(N)

% Plot setting
FontSize_text    = FontSetup.FS_txt;
Fontsize_colorbar = FontSetup.FS_cb;
FN               = FontSetup.FN;

% Coordinates
t_plot = tdata(plot_time_idx) * tr * 60;   % min
x_RV   = (2:N) * Lseg * Lr;
x_P    = (1:N) * Lseg * Lr;

x_RV   = x_RV - mean(x_RV)+Lseg * Lr/2;
x_P    = x_P - mean(x_P)+Lseg * Lr/2;

% =========================
% 1. Intracellular hydraulic pressure
% =========================
figure;
imagesc(x_P, t_plot, P_all(plot_time_idx,:));
set(gca,'YDir','normal');
colormap(cmap);
cb = colorbar;
cb.FontSize = Fontsize_colorbar;
xlabel('AP Position (\mum)');
ylabel('Time (min)');
title('Pressure (kPa)','FontWeight','Normal');
set(gca,'FontSize',FontSize_text,'FontName',FN,'LineWidth',1);

% =========================
% 2. Intracellular osmotic pressure
% =========================
figure;
imagesc(x_RV, t_plot, Pi_all(plot_time_idx,2:N) * Pr);
set(gca,'YDir','normal');
colormap(cmap);
cb = colorbar;
cb.FontSize = Fontsize_colorbar;
xlabel('AP Position (\mum)');
ylabel('Time (min)');
title('Osmolarity (kPa)','FontWeight','Normal');
set(gca,'FontSize',FontSize_text,'FontName',FN,'LineWidth',1);

% =========================
% 5. Interstitial hydraulic ressure
% =========================
figure;
imagesc(x_RV, t_plot, PE_all(plot_time_idx,2:N) * Pr);
set(gca,'YDir','normal');
colormap(cmap);
cb = colorbar;
cb.FontSize = Fontsize_colorbar;
xlabel('AP Position (\mum)');
ylabel('Time (min)');
title('Interstitial Pressure (kPa)','FontWeight','Normal');
set(gca,'FontSize',FontSize_text,'FontName',FN,'LineWidth',1);

end