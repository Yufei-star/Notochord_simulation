function T=tension(R,Rt,S0,k)
h=abs(R)-sqrt(R.^2-Rt^2); %in 3D
S=pi*(Rt^2+h.^2); 
T=k*(S-S0);
end