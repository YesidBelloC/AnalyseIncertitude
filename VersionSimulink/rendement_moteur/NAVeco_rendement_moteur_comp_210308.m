clc, clear all, close all
%%
 SpeedProfileWLTP;
 Profile = WLTPClasse3b;

 Temps   = Profile(:,1);     % s 
 Vitesse = Profile(:,3)/3.6; % m/s
 Accel   = Profile(:,4);     % m/s2 
 Distance = cumtrapz(Temps,Vitesse); % m

%%

Hauteur  = zeros(size(Distance));
Op       = zeros(size(Distance));   % Pente en dégrés
OpPor    = atan(Op/100);            % Pente en %

subplot(3,1,1)
plot(Temps,Distance);
grid on
title('Distance');
xlabel('Temps [s]');
ylabel('Distance [m]');

subplot(3,1,2)
plot(Temps,Vitesse);
grid on
title('Vitesse');
xlabel('Temps [s]');
ylabel('V [m/s]');

subplot(3,1,3)
plot(Temps,Accel);
grid on
title('Acceleration');
xlabel('Temps [s]');
ylabel('Acceleration [m/s2]');
axis([0 max(Temps) -5 6])


%% Calcule des couples et fréquences de rotation
NAVeco_param_EV_210302;
TorqW=[]; rpmW=[];

%ForceHist
FaHist=[];
FwHist=[];
FrHist=[];

% Au niveau de la roue
for i=1:max(size(Temps,1),size(Temps,2))
    TorqW=[TorqW   Rw*(m*Accel(i) + m*g*OpPor(i) + 0.5*rho*Cx*S*Vitesse(i)^2 + m*g*Cr)];

    FaHist=[FaHist 0.5*rho*Cx*S*Vitesse(i)^2];
    FwHist=[FwHist m*g*OpPor(i)];
    FrHist=[FrHist m*g*Cr];
    
    rpmW = [rpmW (Vitesse(i)*30)/(pi*Rw)]; 
end
Ptrac = TorqW/Rw;
Etrac = cumtrapz(Temps,Ptrac)/3600000;

% Au niveau du moteur
TorqM = TorqW/(ig*i0*eff_transm);
rpmM = ig*i0*rpmW;

figure
plot(rpmM,TorqM, 'k*');
grid on
title('Couple moteur selon fréquence de rotation');
xlabel('N [tr/min]');
ylabel('Torq [Nm]');

%% Energie 
if size(TorqM,2)>1, TorqM=TorqM';end
if size(rpmM,2)>1, rpmM=rpmM';end

effm_vec=[]; effr_vec=[];

%**************************************************************************

% Rendements motorisation et régénération variables
effm = 0.9 - Krpm*(rpmM - RPMopt).^2 - Ktorq*(abs(TorqM) - TorqOpt).^2;
effr = 0.74 - Krpm*(rpmM - RPMopt).^2 - Ktorq*(abs(TorqM) - TorqOpt).^2;

%**************************************************************************

Protor = TorqM.*((pi*rpmM)/30);

km=0; kr=0;
for i=1:size(Protor,1)
    if Protor(i)>=0
        km=km+1;
        effm_vec = [effm_vec effm(i)];
        
        eff(i)=1/effm(i);
        eff_vec(i)=effm(i);
    else
        kr=kr+1;
        effr_vec = [effr_vec effr(i)];
        
        eff(i)=effr(i);
        eff_vec(i)=effr(i);
    end
end

Pbat = Protor.*eff';
Ebat = cumtrapz(Temps,Pbat)/3600000;

figure
plot(Temps,Pbat, 'b', 'Linewidth', 1);
title('Puissance électrique à la batterie');
grid on
xlabel('Time [s]');ylabel('Pbat [W]');

figure
plot(Temps,Vitesse*3.6, 'b', 'Linewidth', 1);
grid on
title('Vitesse du véhicule suivant le profil WLTP');
xlabel('Temps [s]');
ylabel('Vitesse [m/s]');
%xlim([0 100]);

figure
plot(Temps,TorqM,'b', 'Linewidth', 1);
grid on
title("Couple appliqué sur l'arbre du rotor");
xlabel('Temps [s]');
ylabel('Couple [Nm]');
%xlim([0 100]);

figure
plot(Temps, Ebat, 'r', 'Linewidth', 1);
title('Energie consommée nette');
grid on
xlabel('Time [s]');
ylabel('Energie [kWh]');
%xlim([0 100]);

figure
plot(Temps,eff_vec, 'r', 'Linewidth', 1);
grid on
title('Rendement de la chaîne de traction');
xlabel('Temps [s]');
ylabel('rendement');
ylim([0 1]);
%xlim([0 100]);

%Energie par Kilometre:
EnergiePer100Km = (Ebat(end)/Distance(end))*100000; % E[kWh]/D[100Km]

disp("Energie consommée kWh       : "+Ebat(end));
disp("Energie consommée kWh/100km : "+EnergiePer100Km);
disp(" ");
disp("rendement moyen de traction   :"+mean(effm_vec));
disp("rendement moyen de régénération :"+mean(effr_vec));