function  [V,h] = all_cell_volume(R,Rt,Lseg)
N=length(R)-1;%# of cells
il=2:N+1;
ip=1:N;
h=(abs(R)-sqrt(R.^2-Rt^2)).*sign(R);
V_SphCap=1/6*pi*h.*(3*Rt^2+h.^2);
V=pi*Rt^2*Lseg+V_SphCap(il)-V_SphCap(ip);
end