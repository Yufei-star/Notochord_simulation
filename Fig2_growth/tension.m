function T=tension(R,Rt,S0,k)
% T=T0*linspace(1,5,length(R));%ones(size(R));

% T=T0*ones(size(R));

% S=2*abs(R).*asin(Rt./abs(R)); %in 2D

h=abs(R)-sqrt(R.^2-Rt^2); %in 3D
S=pi*(Rt^2+h.^2); 
T=k*(S-S0);
end