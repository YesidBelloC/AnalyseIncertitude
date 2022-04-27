function [rMot,rReg] = calcRendement(v,a,m,alpha,S,Cx,Crr,g,rho)

% Vitesse max du véhicule
Vmax_EV = 150/3.6;  % [m/s] 

% Pneu de série 205/45/R17
% Dimensions en m
Haut_flanc = 0.092;
Diam_inter = 0.432;

%Rayon de la roue
Rw = (Diam_inter/2)+Haut_flanc;

% Motorisation
Torq_max = 260;      %[Nm]    à 300 tr/min
RPMmax   = 3673;     %[tr/min]

TorqOpt = (2/5)*Torq_max;
RPMopt  = RPMmax/3;

% Gear ratio
ig = (pi*RPMmax*Rw)/(30*Vmax_EV);   % Rapport de réduction pour Vmax
i0 = 1;                             % rapport de réduction final / Monoréducteur : i0=1
eff_transm = 1;

% Au niveau de la roue
TorqW = Rw.*(m.*a+m.*g.*sin(alpha)+0.5.*rho.*S.*Cx.*v.^2+Crr.*m.*g.*cos(alpha));
rpmW = v.*30./(pi.*Rw);

% Au niveau du moteur
TorqM = TorqW/(ig*i0*eff_transm);
rpmM = ig*i0*rpmW;

% Constantes de la fonction paraboloide du rendement
Krpm  = (1/3)*1e-7;
Ktorq = 1e-5;


% Rendements
rMot = 0.9 - Krpm*(rpmM - RPMopt).^2 - Ktorq*(abs(TorqM) - TorqOpt).^2;
rReg = 0.74 - Krpm*(rpmM - RPMopt).^2 - Ktorq*(abs(TorqM) - TorqOpt).^2;

end