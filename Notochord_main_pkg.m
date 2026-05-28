% -------------------------------------------------------------------------
% Notochord_main_pkg
% Core simulation for perturbation dynamics (e.g., laser ablation)
% -------------------------------------------------------------------------
function [tdata, Geometry, Pressure, Flux] = ...
    Notochord_main_pkg(N, Lseg, alpha_isf, alpha_isf_perturbed, gamma1_perturbed, ...
    alphaw_ext, NumPerturbed, alpha_cell, k, Rt, R, Pr_in, ...
    time, timeperturb_hour, CaseName)

%% ------------------------------------------------------------------------
% Model parameters and physical setup
% ------------------------------------------------------------------------

% Total notochord length (uniform segmentation assumed)
L = N * Lseg;

% Reference membrane area
S0 = pi * Rt^2;

% Store initial geometry
R_init = R;

% -------------------- Transport coefficients ----------------------------
alpha_cell(1) = 0;      % No permeability at left boundary

gamma1 = zeros(1, N);   % Passive ion transport (cell–ECM)
gamma2 = zeros(1, N);   % Active ion transport (cell–ECM)
gamma3 = zeros(1, N);   % Passive ion transport (cell–cell)
gamma4 = zeros(1, N);   % Active ion transport (cell–cell)

% -------------------- Interstitial fluid properties ---------------------
VE   = 5e-2 * ones(1, N);   % Interstitial volume (μm^3)
beta = 100;                 % Pressure–volume coupling
VE0  = 4e-2;                % Preferred interstitial volume

PE  = beta * (VE - VE0);    % Interstitial hydraulic pressure (kPa)
PiE = PE * 800;             % Interstitial osmotic pressure

%% ------------------------------------------------------------------------
% Numerical setup (adaptive time stepping around perturbation)
% ------------------------------------------------------------------------

% Spatial indexing
il = 2:N+1;
ip = 1:N;

% Time discretization
dt_long  = 1e-3;
dt_short = 5e-6;

all_time_points = [ ...
    0:dt_long:(timeperturb_hour - dt_long), ...
    timeperturb_hour:dt_short:(timeperturb_hour + 25*dt_short), ...
    (timeperturb_hour + 26*dt_short):dt_short:time];

numtimestep = length(all_time_points);
tdata = all_time_points;

% Boundary condition:
% 1 → head curvature fixed; 0 → tail curvature fixed
BC_type = 1;

%% ------------------------------------------------------------------------
% Initial state
% ------------------------------------------------------------------------

[V, h] = all_cell_volume(R, Rt, Lseg);
T = tension(R, Rt, S0, k);

% Intracellular hydraulic pressure
P = PE(end) + Pr_in + cumsum(fliplr(2*T(2:end) ./ R(2:end)));
P = fliplr(P);

% Intracellular osmotic pressure
Pi = P + PiE - PE;

% Water flux
J = zeros(1, N);
J(2:N) = alpha_cell(2:N) .* ...
    (P(1:N-1) - P(2:N) - (Pi(1:N-1) - Pi(2:N)));

J_loss = alpha_isf .* (P - PE - Pi + PiE);
J_loss_to_ECM = 0 * alpha_isf;

% Solute flux
Js = -gamma1 .* (PiE - Pi) + gamma2 .* Pi;

Js_cell = zeros(1, N);
Js_cell(2:N) = gamma3(2:N) .* (Pi(1:N-1) - Pi(2:N)) + ...
    gamma4(2:N) .* Pi(1:N-1);

%% ------------------------------------------------------------------------
% Preallocation
% ------------------------------------------------------------------------
R_all      = zeros(numtimestep, N+1);
J_all      = zeros(numtimestep, N);
J_loss_all = zeros(numtimestep, N);
Js_all     = zeros(numtimestep, N);

V_all   = zeros(numtimestep, N);
P_all   = zeros(numtimestep, N);
Pi_all  = zeros(numtimestep, N);

VE_all  = zeros(numtimestep, N);
PE_all  = zeros(numtimestep, N);
PiE_all = zeros(numtimestep, N);

T_all = zeros(numtimestep, N+1);

%% ------------------------------------------------------------------------
% Main simulation loop
% ------------------------------------------------------------------------
for i = 1:numtimestep

    % Adaptive timestep
    if i > 1
        dt = all_time_points(i) - all_time_points(i-1);
    else
        dt = dt_long;
    end

    % Progress display
    disp(strcat(CaseName, ' | progress: ', ...
        num2str(i/numtimestep*100), '% | time: ', ...
        num2str(all_time_points(i)), ' / ', num2str(time)));

    % -------------------- Store current state ----------------------------
    R_all(i,:)      = R;
    J_all(i,:)      = J;
    J_loss_all(i,:) = J_loss;
    Js_all(i,:)     = Js;

    V_all(i,:)  = V;
    P_all(i,:)  = P;
    Pi_all(i,:) = Pi;

    VE_all(i,:)  = VE;
    PE_all(i,:)  = PE;
    PiE_all(i,:) = PiE;

    T_all(i,:) = T;

    % -------------------- Volume dynamics --------------------------------
    S = pi * (Rt^2 + h.^2);

    dV = (J .* S(1:N) ...
        - [J(2:N),0] .* [S(2:N),0] ...
        - J_loss .* [Lseg*(2*pi*Rt)*ones(1,N-1), Lseg*(2*pi*Rt)+S(N+1)]) * dt;

    V = V + dV;

    % Interstitial volume dynamics
    VE = VE + (J_loss .* [Lseg*(2*pi*Rt)*ones(1,N-1), Lseg*(2*pi*Rt)+S(N+1)] ...
        - J_loss_to_ECM .* [Lseg*(2*pi*Rt)*ones(1,N-1), Lseg*(2*pi*Rt)+S(N+1)]) * dt;

    % Darcy flow in interstitial fluid
    VE(2:N)   = VE(2:N)   + alphaw_ext * (PE(1:N-1) - PE(2:N)) * dt;
    VE(1:N-1) = VE(1:N-1) - alphaw_ext * (PE(1:N-1) - PE(2:N)) * dt;

    % -------------------- Geometry update (nonlinear solve) --------------
    h_bc = -(abs(R_init(1)) - sqrt(abs(R_init(1))^2 - Rt^2)) * BC_type ...
        + (abs(R_init(end)) - sqrt(abs(R_init(end))^2 - Rt^2)) * ~BC_type;

    fun = @(h) volume_radius_relation(h, h_bc, Lseg, Rt, V, N, ip, il, BC_type);
    options = optimoptions('fsolve','Display','off');

    if BC_type
        h = [h_bc, fsolve(fun, h(2:end), options)];
    else
        h = [fsolve(fun, h(1:end-1), options), h_bc];
    end

    % Update radius
    R = h/2 + Rt^2 ./ (2*h);

    % -------------------- Pressure and osmolarity ------------------------
    P = PE(end) + Pr_in + cumsum(fliplr(2*T(2:end) ./ R(2:end)));
    P = fliplr(P);

    Pi = ( ...
        -(Js .* [Lseg*(2*pi*Rt)*ones(1,N-1), Lseg*(2*pi*Rt)+S(N+1)]) * dt ...
        + Js_cell .* S(1:N) * dt ...
        - [Js_cell(2:N),0] .* [S(2:N),0] * dt ...
        + Pi_all(i,:) .* V_all(i,:) ...
        ) ./ V;

    PE = beta * (VE - VE0);

    % Fast interstitial ion equilibration assumption
    PiE = PiE_all(1,:);

    % -------------------- Full ion leakage (optional model) --------------
    if gamma1_perturbed == inf
        if all_time_points(i) >= timeperturb_hour
            Pi(NumPerturbed) = ...
                (sum(PiE .* VE) + sum(Pi(NumPerturbed) .* V(NumPerturbed))) ...
                / (sum(VE) + sum(V(NumPerturbed)));

            PiE = ones(size(VE)) * Pi(NumPerturbed(1));
        end
    end

    % -------------------- Update fluxes ----------------------------------
    T = tension(R, Rt, S0, k);

    J = zeros(1, N);
    J(2:N) = alpha_cell(2:N) .* ...
        (P(1:N-1) - P(2:N) - (Pi(1:N-1) - Pi(2:N)));

    J_loss = alpha_isf .* ((P - PE) - (Pi - PiE));
    J_loss_to_ECM = 0 * alpha_isf .* (PE - PE_all(1,1));

    Js = -gamma1 .* (PiE - Pi) + gamma2 .* Pi;

    Js_cell = zeros(1, N);
    Js_cell(2:N) = gamma3(2:N) .* (Pi(1:N-1) - Pi(2:N)) + ...
        gamma4(2:N) .* Pi(1:N-1);

    % -------------------- Perturbation (laser ablation) ------------------
    if all_time_points(i) >= timeperturb_hour

        % Enhanced leakage at perturbed cells
        J_loss(NumPerturbed) = alpha_isf_perturbed * ...
            (P(NumPerturbed) - PE(NumPerturbed));

        % Ion leakage (finite rate)
        if gamma1_perturbed ~= inf
            Js(NumPerturbed) = gamma1_perturbed * ...
                (Pi(NumPerturbed) - PiE(NumPerturbed));
        end

        % Rapid osmotic equilibration
        Pi(NumPerturbed) = PiE(NumPerturbed);

    end

end

%% ------------------------------------------------------------------------
% Output
% ------------------------------------------------------------------------
Geometry = {V_all, R_all, T_all};
Pressure = {P_all, Pi_all, PE_all, PiE_all};
Flux     = {J_all, J_loss_all, Js_all};

end