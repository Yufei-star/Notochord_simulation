function []=plot_config_only(Lseg,N,R,Rt,P,XLIM,Fontsetup)

FN = Fontsetup.FN;
FontSize_colorbar = Fontsetup.FS_cb;% 18;
FontSize_text = Fontsetup.FS_txt;%16;
cmap = viridis(N);

h=(abs(R)-sqrt(R.^2-Rt^2)).*sign(R);
%h(1)=R(1)+sqrt(R(1).^2-Rt^2);
L=sum(Lseg);
numpts=500;

%intracellular pressure and fluid flow
plot([0,L],[2*Rt,2*Rt],'color',[0.28 0.62 0.72],'linewidth',3);hold on;
plot([0,L],[0,0],'color',[0.28 0.62 0.72],'linewidth',3);
theta_span=abs(asin(Rt./R));
for i=1:N+1
    center=sum(Lseg(1:i-1))-R(i)+h(i);
    x_temp=center+R(i)*cos(linspace(-theta_span(i),theta_span(i),numpts));
    y_temp=Rt+R(i)*sin(linspace(-theta_span(i),theta_span(i),numpts));
    x_temp_all(i,:)=x_temp;
    y_temp_all(i,:)=y_temp;
    if i>=2 && ~isempty(P)
        if R(i)*R(i-1)<0
            x=[x_temp_all(i-1,:)';x_temp_all(i,:)'];
            y=[y_temp_all(i-1,:)';y_temp_all(i,:)'];
        else
            x=[x_temp_all(i-1,end:-1:1)';x_temp_all(i,:)'];
            y=[y_temp_all(i-1,end:-1:1)';y_temp_all(i,:)'];
        end
        patch(x,y,P(i-1),'EdgeColor',[1.00 0.78 0.46]);
        hold on;
    end
    %plot(x_temp,y_temp,'k','linewidth',1.5);hold on;
end

colormap(cmap);
c = colorbar;
% clim(CLIM(1,:));
set(c,'FontSize',FontSize_colorbar,'Fontname', FN,'linewidth',1);
set(c,'Location','southoutside');
axis equal
axis off
set(gca,'FontSize',FontSize_text,'Fontname', FN,'linewidth',1);
xlim(XLIM)
ylim([0,2.5*Rt]);
hold off;

end