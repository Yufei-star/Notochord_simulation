function []=plot_config_only(L,N,R,Rt,P,FontSetup)

%Plot setting
FontSize_text = FontSetup.FS_txt;
FontSize_colorbar = FontSetup.FS_cb;
FN = FontSetup.FN;

cmap = viridis(N);

h=(abs(R)-sqrt(R.^2-Rt^2)).*sign(R);
Lseg=L/N;
numpts=500;

%intracellular pressure and fluid flow
plot([0,L],[2*Rt,2*Rt],'color',[0.28 0.62 0.72],'linewidth',3);hold on;
plot([0,L],[0,0],'color',[0.28 0.62 0.72],'linewidth',3);
theta_span=abs(asin(Rt./R));
for i=1:N+1
    center=(i-1)*Lseg-R(i)+h(i);
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
end

colormap(cmap);
c = colorbar;
set(c,'FontSize',FontSize_colorbar,'Fontname', FN,'linewidth',1);
set(c,'Location','southoutside');
axis equal
axis off
set(gca,'FontSize',FontSize_text,'Fontname', FN,'linewidth',1);
xlim([0-Lseg L+Lseg])
ylim([0,2.5*Rt]);
hold off;

end