function []=Plot1_dynamical_curve_representation_Intra(tdata,R_all,V_all,T_all,...
    P_all,Pi_all,PE_all,PiE_all,...
    J_all,J_loss_all,Js_all,...
    Pr,tr,Lr,Lseg,FontSetup)
% With unit, dynamical, curve form same as what Parsa plotted

N = length(V_all(1,:));

%Plot setting
FontSize_text = FontSetup.FS_txt;
Fontsize_colorbar = FontSetup.FS_cb;
FN = FontSetup.FN;
cbstring='Position (\times10^2μm)';

Cneg=[255,0,255]; Cmed=[128,128,128]; Cpos=[0,255,255];
cmap=[[linspace(Cneg(1),Cmed(1),N/2),linspace(Cmed(1),Cpos(1),N/2)]',...
    [linspace(Cneg(2),Cmed(2),N/2),linspace(Cmed(2),Cpos(2),N/2)]',...
    [linspace(Cneg(3),Cmed(3),N/2),linspace(Cmed(3),Cpos(3),N/2)]']/255;

cmap = viridis(N);

% Plot Volume dynamics (percentage)
figure;
for i=1:N-1
    plot(tdata*tr*60,100*(V_all(:,i+1)'-V_all(1,i+1))/V_all(1,i+1),'Color',cmap(i+1,:),'LineWidth', 2);
    hold on;
end
colormap(cmap);
cb = colorbar;
caxis([-(N/2-1)*Lseg*Lr, N/2*Lseg*Lr]/1e2);
xlim([min(tdata),max(tdata)]*tr*60);
ylim([-60,5]);
xlabel('Time (min)');
ylabel('Volume Change (%)');
set(gca,'FontSize',FontSize_text,'Fontname', FN,'linewidth',1);
cb = colorbar;
cb.Label.String = cbstring;
cb.Label.Position = [2.5 0 0];  % [x y z] relative to colorbar center
cb.Label.Rotation = 90;         % Keep it vertical (default for right-side bar)
cb.FontSize = Fontsize_colorbar;

% Plot curvature dynamics
figure;
for i=1:N-1
    plot(tdata*tr*60,(1./R_all(:,i+1)'-1/R_all(1,i+1))/Lr*1e3,'Color',cmap(i+1,:),'LineWidth', 2);
    hold on;
end
colormap(cmap);
cb = colorbar;
cb.Label.String = cbstring;
cb.Label.Position = [2.5 0 0];  % [x y z] relative to colorbar center
cb.Label.Rotation = 90;         % Keep it vertical (default for right-side bar)
cb.FontSize = Fontsize_colorbar;
caxis([-(N/2-1)*Lseg*Lr, N/2*Lseg*Lr]/1e2);
xlim([min(tdata),max(tdata)]*tr*60);
ylim([-0.082,0.02]*1e3)
xlabel('Time (min)');
ylabel('\DeltaC (\times10^{-3}μm^{-1})');
set(gca,'FontSize',FontSize_text,'Fontname', FN,'linewidth',1);

% Plot tension dynamics
figure;
for i=1:N
    plot(tdata*tr*60,T_all(:,i)'/T_all(1,i),'Color',cmap(i,:),'LineWidth', 2);
    hold on;
end
d_values = linspace(-N/2*Lseg*Lr, N/2*Lseg*Lr, N);  % different distance to punture site values
colormap(cmap);
cb = colorbar;
cb.Label.String = cbstring;
cb.Label.Position = [2.5 0 0];  % [x y z] relative to colorbar center
cb.Label.Rotation = 90;         % Keep it vertical (default for right-side bar)
cb.FontSize = Fontsize_colorbar;
caxis([-N/2*Lseg*Lr, N/2*Lseg*Lr]/1e2);
xlim([min(tdata),max(tdata)]*tr*60);
xlabel('Time (min)');
ylabel('Tension Fold-change');
set(gca,'FontSize',FontSize_text,'Fontname', FN,'linewidth',1);

end