function []=plot_notochord_video(Lseg,N,R,Rt,t,Q,Q1,Q2,Q3,J,tr,CLIM,CbOri,XLIM,Fontsetup)

FN = Fontsetup.FN;
FontSize_colorbar = Fontsetup.FS_cb;% 12;
FontSize_text = Fontsetup.FS_txt;%16;

cmap = viridis(N);
YLIM = [-3*Rt,3*Rt];

h=(abs(R)-sqrt(R.^2-Rt^2)).*sign(R);
L=sum(Lseg);
numpts=500;

%intracellular pressure and fluid flow
subplot(2,1,1);
plot([0,L],[Rt,Rt],'color',[0.28 0.62 0.72],'linewidth',3);hold on;
plot([0,L],[-Rt,-Rt],'color',[0.28 0.62 0.72],'linewidth',3);
theta_span=abs(asin(Rt./R));
for i=1:N+1
    center=sum(Lseg(1:i-1))-R(i)+h(i);
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
quiver(linspace(Lseg(1)/2,L-Lseg(N)/2,N),Rt,zeros(1,N),J,'r','linewidth',2,'AutoScale', 'off');

colormap(cmap);
c = colorbar;
clim(CLIM(1,:));
set(c,'FontSize',FontSize_colorbar,'Fontname', FN,'linewidth',1);
if ~CbOri
    set(c,'Location','southoutside');
end
axis equal
axis off
title('$P-P^{(e)}$ (kPa)','FontWeight', 'Normal','Interpreter','latex');
set(gca,'FontSize',FontSize_text,'Fontname', FN,'linewidth',1);
xlim(XLIM);
ylim(YLIM);
hold off;

%intracellular osmotic pressure
subplot(2,1,2);
plot([0,L],[Rt,Rt],'color',[0.28 0.62 0.72],'linewidth',3);hold on;
plot([0,L],[-Rt,-Rt],'color',[0.28 0.62 0.72],'linewidth',3);
theta_span=abs(asin(Rt./R));
for i=1:N+1
    center=sum(Lseg(1:i-1))-R(i)+h(i);
    x_temp=center+R(i)*cos(linspace(-theta_span(i),theta_span(i),numpts));
    y_temp=Rt+R(i)*sin(linspace(-theta_span(i),theta_span(i),numpts));
    x_temp_all(i,:)=x_temp;
    y_temp_all(i,:)=y_temp;
    if i>=2 && ~isempty(Q1)
        if R(i)*R(i-1)<0
            x=[x_temp_all(i-1,:)';x_temp_all(i,:)'];
            y=[y_temp_all(i-1,:)';y_temp_all(i,:)'];
        else
            x=[x_temp_all(i-1,end:-1:1)';x_temp_all(i,:)'];
            y=[y_temp_all(i-1,end:-1:1)';y_temp_all(i,:)'];
        end
        patch(x,y-Rt,Q1(i-1),'EdgeColor',[1.00 0.78 0.46]);
        hold on;
    end
end

colormap(cmap);
c = colorbar;
clim(CLIM(2,:));
set(c,'FontSize',FontSize_colorbar,'Fontname', FN,'linewidth',1);
axis equal
axis off
title('$\Pi-\Pi^{(e)}$ (kPa)','FontWeight', 'Normal','Interpreter','latex');
set(gca,'FontSize',FontSize_text,'Fontname', FN,'linewidth',1);
xlim(XLIM);
ylim(YLIM);
hold off;

sgtitle(strcat('Time: ',num2str(t*tr),' h'),'FontSize',20);

end