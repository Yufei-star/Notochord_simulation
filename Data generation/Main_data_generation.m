%% PART 1,2: Base value (Fig3 and Sup3)
clc;clear;
get_notochord_params_base;
CaseName='BaseValueResult';
[tdata,Geometry,Pressure,Flux]=Notochord_main_pkg(N,Lseg,alpha_isf,alpha_isf_perturbed,gamma1_perturbed,...
                                                  alphaw_ext,NumPerturbed,alpha_cell,k,Rt,R,Pr_in,...
                                                  time,timeperturb_hour,CaseName); %tdata: hour
save(CaseName);

%% PART 3-1: Assymetrical (Fig 4)
% High curvature case is the same data as Fig 3 (BaseValueResult)
clc;clear;
get_notochord_params_base;

% Initial configuration
R0 = 2 * Rt;%1.5; % Initial radius
% configuration 1 
R = R0 * ones(1, N+1); % Boundary radii
R(1) = -R0; % Left boundary
kp = 1 ./ R;

CaseName='LowCurv_R=1';
[tdata,Geometry,Pressure,Flux]=Notochord_main_pkg(N,Lseg,alpha_isf,alpha_isf_perturbed,gamma1_perturbed,...
                                                  alphaw_ext,NumPerturbed,alpha_cell,k,Rt,R,Pr_in,...
                                                  time,timeperturb_hour,CaseName); %tdata: hour
save(CaseName);

%% PART 3-2: Initial mean curv vs. response asymmetry (Fig 5A)
% Parameter scan for generating and saving data
clc;clear;
get_notochord_params_base;
C0_all=[linspace(-1,-0.0667,11),linspace(0.0667,1,11)]; %Initial curvature
R0_base = 1.25*Rt;
R0_all = 1./C0_all;

%Strating parameter scan
for i = 1:length(R0_all)
    CaseName=strcat('RI0_',num2str(R0_all(i)));
    R_temp = R0_base * ones(1, N+1); % Boundary radii
    R_temp(1) = -R0_base; % Left boundary
    R_temp(NumPerturbed(1)-10:NumPerturbed(end)+10) = R0_all(i);
    [tdata,Geometry,Pressure,Flux]=Notochord_main_pkg(N,Lseg,alpha_isf,alpha_isf_perturbed,gamma1_perturbed,...
        alphaw_ext,NumPerturbed,alpha_cell,k,Rt,R_temp,Pr_in,time,timeperturb_hour,CaseName); %tdata: hour
    save(strcat('C0_',num2str(C0_all(i)),'.mat'));
end
%% PART 3-3: Initial curv asymmetry vs. response asymmetry (Fig 5B)
% Parameter scan for generating and saving data (local curvature variation)
clc;clear;
get_notochord_params_base;
% Initial configuration
Rt=0.5;
R0_base=1.25*Rt;
C0_left_all = linspace(1.05,1.95,12);
C0_right_all = 3 - C0_left_all; % Ensure the same mean curvature while allowing curvature assymetry
Left_idx = NumPerturbed(1)-10:NumPerturbed(1); %local curvature varations in left 10 cells to puncture site
Right_idx = NumPerturbed(end):NumPerturbed(end)+10; %local curvature varations in right 10 cells to puncture site

%Strating parameter scan
% All posterior-facing
for i = 1:length(C0_left_all)
    CaseName=strcat('CI_left0_',num2str(C0_left_all(i)));
    R_temp = R0_base * ones(1, N+1); % Boundary radii\
    R_temp(1) = -R0_base; % Left boundary
    R_temp(Left_idx) = 1/C0_left_all(i); % Left 10 cells to puncture
    R_temp(Right_idx) = 1/C0_right_all(i); % Right 10 cells to puncture
    [tdata,Geometry,Pressure,Flux]=Notochord_main_pkg(N,Lseg,alpha_isf,alpha_isf_perturbed,gamma1_perturbed,...
        alphaw_ext,NumPerturbed,alpha_cell,k,Rt,R_temp,Pr_in,time,timeperturb_hour,CaseName); %tdata: hour
    save(strcat('C0_left_',num2str(C0_left_all(i)),'.mat'));
end

% Opposing curvature profile
for i = 1:length(C0_left_all)
    CaseName=strcat('CI_left0_',num2str(C0_left_all(i)));
    R_temp = R0_base * ones(1, N+1); % Boundary radii\
    R_temp(1) = -R0_base; % Left boundary
    R_temp(Left_idx) = 1/C0_left_all(i); % Left 10 cells to puncture
    R_temp(Right_idx) = -1/C0_right_all(i); % Right 10 cells to puncture
    [tdata,Geometry,Pressure,Flux]=Notochord_main_pkg(N,Lseg,alpha_isf,alpha_isf_perturbed,gamma1_perturbed,...
        alphaw_ext,NumPerturbed,alpha_cell,k,Rt,R_temp,Pr_in,time,timeperturb_hour,CaseName); %tdata: hour
    save(strcat('Neg_C0_left_',num2str(C0_left_all(i)),'.mat'));
end

% All anterior-facing
for i = 1:length(C0_left_all)
    CaseName=strcat('CI_left0_',num2str(C0_left_all(i)));
    R_temp = R0_base * ones(1, N+1); % Boundary radii\
    R_temp(1) = -R0_base; % Left boundary
    R_temp(Left_idx) = -1/C0_left_all(i); % Left 10 cells to puncture
    R_temp(Right_idx) = -1/C0_right_all(i); % Right 10 cells to puncture
    [tdata,Geometry,Pressure,Flux]=Notochord_main_pkg(N,Lseg,alpha_isf,alpha_isf_perturbed,gamma1_perturbed,...
        alphaw_ext,NumPerturbed,alpha_cell,k,Rt,R_temp,Pr_in,time,timeperturb_hour,CaseName); %tdata: hour
    save(strcat('Neg2_C0_left_',num2str(C0_left_all(i)),'.mat'));
end

%% PART 4-1: Volume change with different membrane elasticity (Fig 3H)
clc;clear;
get_notochord_params_base;
k_all = [0.1 0.25 0.5]; % Varying elasticity
for i=1:length(k_all)
    % Initial configuration
    CaseName=strcat('Vloss_Elasticity_',num2str(k_all(i)));
    [tdata,Geometry,Pressure,Flux]=Notochord_main_pkg(N,Lseg,alpha_isf,alpha_isf_perturbed,gamma1_perturbed,...
        alphaw_ext,NumPerturbed,alpha_cell,k_all(i),Rt,R,Pr_in,...
        time,timeperturb_hour,CaseName); %tdata: hour
    save(CaseName);
end

%% PART 4-2: Volume change with different Pressure differnece (Fig 5E)
% The pressure difference is gnerated by varying mean curvature of the
% system (global)
clc;clear;
get_notochord_params_base;
C0_all=linspace(0.0667,1.25,14);
R0_all = 1./C0_all;
%Strating parameter scan
for i = 1:length(R0_all)
    CaseName=strcat('RI0_',num2str(R0_all(i)));
    R_temp = R0_all(i) * ones(1, N+1); % Boundary radii
    R_temp(1) = -R0_all(i); % Left boundary
    [tdata,Geometry,Pressure,Flux]=Notochord_main_pkg(N,Lseg,alpha_isf,alpha_isf_perturbed,gamma1_perturbed,...
        alphaw_ext,NumPerturbed,alpha_cell,k,Rt,R_temp,Pr_in,time,timeperturb_hour,CaseName); %tdata: hour
    save(strcat('C0_',num2str(C0_all(i)),'.mat'));
end

%% PART 5: Decaying profile for different membrane permeability (Fig S5 A-B)
% For different initial configuration, the data is taken from PART3
clc;clear;
get_notochord_params_base;
alpha_isf_all = [0 0.001 0.01 0.1];
for i = 1:length(alpha_isf_all)
    CaseName=strcat('2alpha_isf_',num2str(alpha_isf_all(i)));
    alpha_isf_temp = alpha_isf_all(i)*ones(1, N);
    [tdata,Geometry,Pressure,Flux]=Notochord_main_pkg(N,Lseg,alpha_isf_temp,alpha_isf_perturbed,gamma1_perturbed,...
        alphaw_ext,NumPerturbed,alpha_cell,k,Rt,R,Pr_in,time,timeperturb_hour,CaseName); %tdata: hour
    save(strcat('alpha_isf_',num2str(alpha_isf_all(i)),'.mat'));
end