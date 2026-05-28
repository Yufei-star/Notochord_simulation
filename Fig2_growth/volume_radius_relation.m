function F = volume_radius_relation(h,h_bc,Lseg,Rt,V,N,ip,il,BC_type)
%BC_type=1: head cell memb curvature is fixed; BC_type=2: tail cell memb
%curvature is fixed
h_temp=[h_bc,h]*BC_type+[h,h_bc]*~BC_type;
V_SphCap=1/6*pi*h_temp.*(3*Rt^2+h_temp.^2);
F=zeros(N,1);
F=pi*Rt^2*Lseg+V_SphCap(il)-V_SphCap(ip)-V;
% max(abs(F(:)))
end