function []=plot_notochord_video(L,N,R,Rt,t,Q,Q1,Q2,Q3,J,tr,CLIM,CbOri,FontSetup)

FontSize_colorbar = FontSetup.FS_cb;%12;
FontSize_text = FontSetup.FS_txt;%16;
FN = FontSetup.FN;
cmap = viridis(N);
YLIM = [-3*Rt,3*Rt];

h=(abs(R)-sqrt(R.^2-Rt^2)).*sign(R);
Lseg=L/N;
numpts=500;

%intracellular pressure and fluid flow
subplot(2,1,1);
plot([0,L],[Rt,Rt],'color',[0.28 0.62 0.72],'linewidth',3);hold on;
plot([0,L],[-Rt,-Rt],'color',[0.28 0.62 0.72],'linewidth',3);
theta_span=abs(asin(Rt./R));
for i=1:N+1
    center=(i-1)*Lseg-R(i)+h(i);
    x_temp=center+R(i)*cos(linspace(-theta_span(i),theta_span(i),numpts));
    y_temp=Rt+R(i)*sin(linspace(-theta_span(i),theta_span(i),numpts));
    x_temp_all(i,:)=x_temp;
    y_temp_all(i,:)=y_temp;
    if i>=2 && ~isempty(Q)
        if R(i)*R(i-1)<0
            x=[x_temp_all(i-1,:)';x_temp_all(i,:)'];
            y=[y_temp_all(i-1,:)';y_temp_all(i,:)'];
        else
            x=[x_temp_all(i-1,end:-1:1)';x_temp_all(i,:)'];
            y=[y_temp_all(i-1,end:-1:1)';y_temp_all(i,:)'];
        end
        patch(x,y-Rt,Q(i-1),'EdgeColor',[1.00 0.78 0.46]);
        hold on;
    end
end
quiver(linspace(Lseg/2,L-Lseg/2,N),Rt,zeros(1,N),J,'r','linewidth',2,'AutoScale', 'off');

colormap(cmap);
c = colorbar;
set(c,'FontSize',FontSize_colorbar,'Fontname', FN,'linewidth',1);
if ~CbOri
    set(c,'Location','southoutside');
end
axis equal
axis off
title('Intracellular hydraulic pressure (kPa)','FontWeight', 'Normal');
set(gca,'FontSize',FontSize_text,'Fontname', FN,'linewidth',1);
xlim([0-Lseg L+Lseg])
ylim(YLIM);
hold off;

%interstitial pressure and fluid flow
subplot(2,1,2);
numpts2=4;
x_temp_all=[];
y_temp_all=[];
plot([0,L],[Rt,Rt],'color',[0.28 0.62 0.72],'linewidth',3);hold on;
plot([0,L],[-Rt,-Rt],'color',[0.28 0.62 0.72],'linewidth',3);
[xmesh,ymesh]=meshgrid(linspace(0,L,length(Q2)-1),linspace(0,2*Rt,numpts2));
J_temp=diff(Q2);
v=zeros(size(xmesh));
u=J_temp.*ones(size(xmesh));
for i=1:N+1
    x_temp=(i-1)*Lseg*ones(numpts2,1);
    y_temp=linspace(0,2*Rt,numpts2)';
    x_temp_all(:,i)=x_temp;
    y_temp_all(:,i)=y_temp;
    if i>=2 && ~isempty(Q2)
        x=[x_temp_all(end:-1:1,i-1);x_temp_all(:,i)];
        y=[y_temp_all(end:-1:1,i-1);y_temp_all(:,i)];
        patch(x,y-Rt,Q2(i-1),'EdgeColor','None');
        hold on;
    end
end

colormap(cmap);
c = colorbar;
clim(CLIM(3,:));
set(c,'FontSize',FontSize_colorbar,'Fontname', FN,'linewidth',1);
if ~CbOri
    set(c,'Location','southoutside');
end
axis equal
axis off
title('Interstitial hydraulic pressure (kPa)','FontWeight', 'Normal');
set(gca,'FontSize',FontSize_text,'Fontname', FN,'linewidth',1);
xlim([0-Lseg L+Lseg])
ylim(YLIM);
hold off;

sgtitle(strcat('Time: ',num2str(t*60*tr),' min'));

end