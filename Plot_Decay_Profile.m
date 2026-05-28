function []=Plot_Decay_Profile(tdata,Geometry,Pressure,Flux,timeperturb_hour,NumPerturbed,Lr,Lseg,FigNum,FontSetup)

FontSize_text = FontSetup.FS_txt;
Fontsize_colorbar = FontSetup.FS_cb;
FN = FontSetup.FN;

V_all = Geometry{1}; R_all = Geometry{2}; T_all = Geometry{3};
P_all = Pressure{1}; Pi_all = Pressure{2};
PE_all = Pressure{3}; PiE_all = Pressure{4};
J_all = Flux{1}; J_loss_all = Flux{2}; Js_all = Flux{3};
C_all = 1./R_all;

NumSelect = NumPerturbed(end)-10:NumPerturbed(end)+10;

[~,idx]=min(abs(tdata-timeperturb_hour));
dVm_all=max(abs(V_all(idx:end,:)-V_all(1,:))./V_all(1,:));
dV_mag=dVm_all(NumSelect)/dVm_all(NumPerturbed(1));
dCm_all=max(abs(C_all(idx:end,:)-C_all(1,:)))/Lr;
dC_mag=dCm_all(NumSelect(2:end))/dCm_all(NumPerturbed(end)+1);

xx = NumSelect*Lseg*Lr;

dV_mag_temp = dV_mag;
dV_mag_temp([10 11]) = [];
xx_temp = xx;
xx_temp([10 11]) = [];
figure(FigNum(1));
plot(xx_temp-1220,dV_mag_temp,'linewidth',2);
xlabel('Distance from Puncture (μm)');
ylabel('Relative \DeltaV');
xlim([min(xx) max(xx)]-1220);
ylim([-0.02 1.05]);
set(gca,'FontSize',FontSize_text,'Fontname', FN,'linewidth',1);
hold on;

xx(10)=[];
dC_mag(10)=[];
figure(FigNum(2));
plot(xx(2:end)-1230,dC_mag,'linewidth',2);
hold on;
xlabel('Distance from Puncture (μm)');
ylabel('Relative \DeltaC');
xlim([min(xx(2:end)) max(xx(2:end))]-1230);
ylim([-0.02 1.05]);
set(gca,'FontSize',FontSize_text,'Fontname', FN,'linewidth',1);

end


function [xmin, xmax] = find_decay_xlim(x, y)

% Smooth slightly to avoid noisy derivative triggering false plateau
y_s = smoothdata(y, 'movmean', 5);

% x at maximum y
[~, idx_max] = max(y_s);
xmin = x(idx_max);

% Numerical slope
dy = gradient(y_s)./gradient(x);

% Threshold for "plateau"
thr = 0.02 * max(abs(dy));   % can tune this
n_consecutive = 5;           % require plateau for several points

idx_plateau = [];

for i = idx_max+1 : length(y)-n_consecutive+1
    if all(abs(dy(i:i+n_consecutive-1)) < thr)
        idx_plateau = i;
        break;
    end
end

if isempty(idx_plateau)
    xmax = x(end);
else
    xmax = x(idx_plateau);
end

% Safety: avoid invalid xlim if plateau too close to max
if xmax <= xmin
    xmax = x(end);
end

end