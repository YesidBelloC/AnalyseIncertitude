lambda0 = 10e-3;
lambda2 = 4e-8;
rho     = 1.25;      % Masse volumique de l'air
Cr      = 0.01;      % Coefficient de résistance au roulement
g       = 9.81;      % Accélération gravitationnelle

rend_trac=0.8;        % (p.u)
rend_rege=0.64;       % (p.u)

%% Peugeot e-208 (électrique)
m   = 1500;         % [kg]   Masse
SCx = 0.62;
Cx  = 0.29;         % Coefficient de traînée

S = SCx/Cx;         % [m^2]  Surface frontale
Vmax_EV = 150/3.6;  % [m/s] Vitesse max du véhicule

% Pneu de série 205/45/R17
% Dimensions en m
Larg       = 0.205;
Haut_flanc = 0.092;
Diam_inter = 0.432;

%Rayon de la roue
Rw = (Diam_inter/2)+Haut_flanc;

% Motorisation
Pmax     = 100000;   %[W]     à 3673 tr/min
Torq_max = 260;      %[Nm]    à 300 tr/min
RPMmax   = 3673;     %[tr/min]

TorqOpt = (2/5)*Torq_max;
RPMopt  = RPMmax/3;

% Gear ratio
ig = (pi*RPMmax*Rw)/(30*Vmax_EV);   % Rapport de réduction pour Vmax
i0 = 1;                             % rapport de réduction final / Monoréducteur : i0=1
eff_transm = 1;

% Constantes de la fonction paraboloide du rendement
Krpm  = (1/3)*1e-7;
Ktorq = 1e-5;