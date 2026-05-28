% Base values for all parameters
N = 120;       % Number of cells (Should be larger 100)
Lseg = 0.5; % Length of each cell
k = 0.5;        % Elastic modulus of the membrane (for new version with elastic membrane)
alpha_isf = 1e-2 * ones(1, N); % Permeability of lateral sheath !!!!!!!!!!!!!!!!!!!!!!1 1e-2
alpha_cell = zeros(1, N); % Permeability between cells (left membrane)
alphaw_ext=2; % water transport in interstitial fluid - Darcy's law
alpha_isf_perturbed=100;%100;%5e-1;%alpha_isf(NumPerturbed)*10;%10 !!!!!!!!!!!!!!!!!!2 was 5e-1 before
gamma1_perturbed=1;%inf;
NumPerturbed=[N/2,N/2+1];
Pr_in=0; %reference intracellular pressure
time=0.125;
timeperturb_hour=0.025;% perturbation @ time point (h)
Idx_select_cell=NumPerturbed(1)-10:NumPerturbed(2)+10;
Idx_select_mem = [Idx_select_cell,Idx_select_cell(end)+1];

% Initial configuration
Rt=0.5;
R0 = 1.25 * Rt;%1.5; % Initial radius
% configuration 1 
R = R0 * ones(1, N+1); % Boundary radii
R(1) = -R0; % Left boundary
kp = 1 ./ R;

% Scale parameters
Lr=2*20; %2Rt (um)
Pr=1; % kPa
tr=1*2; % h