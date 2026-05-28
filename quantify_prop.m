% Quantifying spatical decaying and local time scale

function [t_dec,d_dec]=quantify_prop(tdata,V_all,C_all,NumPert,tPert,Lr,tr) %decaying time scale t_dec and spatialy decaying distance d_dec

% Local time scales
[~,idx]=min(abs(tdata-tPert));
V_loc = V_all(idx:end,NumPert(1))-V_all(end,NumPert(1));
C_loc = C_all(idx:end,NumPert(1)+1)-C_all(end,NumPert(1)+1);

% Suppose you have your data vectors:
t = tdata(idx:end)-tdata(idx);  % time data
fv = fit(t', V_loc, 'exp2'); %double exponential: a1*exp(-t/t1)+a2*exp(-t/t2)
fc = fit(t', C_loc, 'exp2'); %double exponential: a1*exp(-t/t1)+a2*exp(-t/t2)
t1v=-1/fv.b; t2v=-1/fv.d; t_temp=t1v;
a1v=fv.a; a2v=fv.c;a_temp=a1v;
if t1v>t2v
t1v=t2v; t2v=t_temp;
a1v=fv.c;a2v=a_temp;
end
disp(['ratio_V:',num2str(a2v/a1v)]);

t1c=-1/fc.b; t2c=-1/fc.d; t_temp=t1c;
a1c=fc.a; a2c=fc.c;a_temp=a1c;
if t1c>t2c
t1c=t2c; t2c=t_temp;
a1c=fc.c;a2c=a_temp;
end
disp(['ratio_C:',num2str(a2c/a1c)]);

t_dec=[t1v,t2v;t1c,t2c]*tr*60;

dVm_all=max(abs(V_all(idx:end,:)-V_all(1,:))./V_all(1,:));
dV_mag=dVm_all/max(dVm_all);
dCm_all=max(abs(C_all(idx:end,:)-C_all(1,:))./C_all(1,:));
dC_mag=dCm_all/max(dCm_all);

dVL_mag=fliplr(dV_mag(1:NumPert(1)));
dVR_mag=dV_mag(NumPert(end):end);
dCL_mag=fliplr(dC_mag(1:NumPert(1)));
dCR_mag=dC_mag(NumPert(end)+1:end);

fVL = fit((1:NumPert(1))'-1, dVL_mag', 'exp1');
fVR = fit((1:length(dVR_mag))'-1, dVR_mag', 'exp1');
fCL = fit((1:NumPert(1))'-1, dCL_mag', 'exp1');
fCR = fit((1:length(dCR_mag))'-1, dCR_mag', 'exp1');

d_decay_VL = -1/fVL.b;
d_decay_VR = -1/fVR.b;
d_decay_CL = -1/fCL.b;
d_decay_CR = -1/fCR.b;

d_dec=[d_decay_VL,d_decay_VR;d_decay_CL,d_decay_CR]*Lr;

%Check the fitting
% plot(fv,t,V_loc);
% set(gca,'FontSize',24,'Fontname', 'Arial','linewidth',1);
% xlabel('Time (h)');
% ylabel('Volume (a.u.)');
% figure;
% plot(fc,t,C_loc);
% xlabel('Time (h)');
% ylabel('Curvature (a.u.)');
% set(gca,'FontSize',24,'Fontname', 'Arial','linewidth',1);
% figure;
% plot(fVL,(1:NumPert(1))-1,dVL_mag);
% xlabel('Position');
% ylabel('Volume (L)');
% set(gca,'FontSize',24,'Fontname', 'Arial','linewidth',1);
% figure;
% plot(fVR,(1:length(dVR_mag))-1, dVR_mag);
% xlabel('Position');
% ylabel('Volume (R)');
% set(gca,'FontSize',24,'Fontname', 'Arial','linewidth',1);
% figure;
% plot(fCL,(1:NumPert(1))-1,dCL_mag);
% xlabel('Position');
% ylabel('Curvature (L)');
% set(gca,'FontSize',24,'Fontname', 'Arial','linewidth',1);
% figure;
% plot(fCR,(1:length(dCR_mag))-1, dCR_mag);
% xlabel('Position');
% ylabel('Curvature (R)');
% set(gca,'FontSize',24,'Fontname', 'Arial','linewidth',1);

end