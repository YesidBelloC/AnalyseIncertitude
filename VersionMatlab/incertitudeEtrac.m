%Permet de calculer les incertitudes des forces, de la puissance et de
%l energie, le poids des parametres qui les influence et l incertitude
%maximale de chaque parametre pour ne pas depasser un pourcentage
%d incertitude sur l energie
%% Initialisation

clear;
close all;
clc;

addpath(genpath('../03_incertitude_quadratique'));

%% Parametres pour le calcul des incertitudes et des poids

%Recuperation temps,distance,vitesse et acceleration du cycle WLTP
[t,d,v,a] = getWLTP();  

%Masse moyenne véhicule en kg
m = 1240;                   %Masse moyenne d'un vehicule

%Pente de la route en %
alpha = 0;                  %Pente nulle

%Coefficients aérodynamiques
S = 2.3;                    %Maitre couple moyen
Cx = 0.36;                  %Coefficient de traine moyen

%Coefficient de resistance au roulement 
lambda0 = 0.005;
lambda1 = 1000;
lambda2 = 1.2312; 
p = 200000;                 %Pression moyenne en Pascal
CrrMod = 0;
CrrNom = 0.01;              %0 : Coefficient de resistance au roulement moyen
CrrCalc = calcCrr(lambda0,lambda1,lambda2,p,v);  %1 : Crr =lambda0+1./p.*(lambda1+lambda2.*v.^2)
if (CrrMod == 0)
    Crr = CrrNom;
else
    Crr = CrrCalc;
end

%Intensite gravitationnelle en m/s²
h = 0;                      %Altitude en metre
gMod = 0;
gNom = 9.81;                %0 : Intensite gravitationnelle moyenne
gCalc = calcG(h);         	%1 : g = (G*Mt)./((Rt+h).^2)
if (gMod == 0)
    g = gNom;
else
    g = gCalc;
end

%Masse volumique de l'aire en kg/m³
T = 285;                    %Température en Kelvin
P = 102900;                 %Pression atmospherique en Pascal
rhoMod = 0;
rhoNom = 1.292;             %0 : Masse volumique moyenne de l'air
rhoHT = calcRhoHT(h,T,g);   %1 : rho = P0.*M.*exp(-h.*M.*g./(R.*T))./(R.*T);
rhoPT = calcRhoPT(P,T);     %2 : rho = P.*M./(R.*T);
if (rhoMod == 0)
    rho = rhoNom;
elseif (rhoMod == 1)
    rho = rhoHT;
else
    rho = rhoPT;
end

%Rendements variables
%rMot:rendement en phase motrice / rReg:rendement en phase de regeneration
%[rMot,rReg] = calcRendement(v,a,m,alpha,S,Cx,Crr,g,rho);
% Rendements pour une voiture electrique
rMot = 0.8216.*ones(size(t));
rReg = 0.6774.*ones(size(t));

%Pourcentage d'incertitude maximum sur Etrac pour le calcul des
%incertitudes maximales
pE = 1;                    %(%)

%Incertitude des paramètres
uh = 0;            %(m)
uT = 0;            %(K)
uP = 0;                    %(P)
uSCx = 0; 
um = 10;           %(kg)
up = 0;         %(P)
ualpha = 0;          %(rad)
uv = 0;              %(m/s)

%% Calculs energie, puissance, incertitudes et poids

%Energie (kWh) et puissance (W) de traction / Incertitudes forces, puissance et energie de traction
[uFaero,pFaero,uFrr,pFrr,uFgrav,pFgrav,uFnet,pFnet,uPtrac,pPtrac,uEtrac,uMax,Ptrac,Etrac] = calcEtracIncertitudesPoids(t,v,a,m,alpha,S,Cx,lambda0,lambda1,lambda2,p,Crr,h,g,T,P,rhoMod,rho,uT,uP,uh,um,up,uv,ualpha,uSCx,pE,rMot,rReg);

%Poids des parametres de Faero
prhoFaero = pFaero(1);pSCxFaero = pFaero(2);pvFaero = pFaero(3);
%Poids des parametres de Frr
pCrrFrr = pFrr(1);pmFrr = pFrr(2);pgFrr = pFrr(3);palphaFrr = pFrr(4);
%Poids des paramètres de Fgrav
pmFgrav = pFgrav(1);pgFgrav = pFgrav(2);palphaFgrav = pFgrav(3);
%Poids des parametres de Fnet
pmFnet = pFnet(1);
%Poids des parametres de la puissance de traction
phP = pPtrac(1);pTP = pPtrac(2);pPP = pPtrac(3);pSCxP = pPtrac(4);pmP = pPtrac(5);ppP = pPtrac(6);palphaP = pPtrac(7);pvP = pPtrac(8);
%Incertitudes maximales des paramètres
uhMax = uMax(1);uTMax = uMax(2);uPMax = uMax(3);uSCxMax = uMax(4);umMax = uMax(5);upMax = uMax(6);ualphaMax = uMax(7);uvMax = uMax(8);
%Incertitude moyenne des forces
uFaeroMoy = mean(uFaero);uFrrMoy = mean(uFrr);uFgravMoy = mean(uFgrav);uFnetMoy = mean(uFnet);uFtracMoy = uFaeroMoy+uFrrMoy+uFgravMoy+uFnetMoy;

PtracMoy = mean(Ptrac);
uPtracMoy = mean(uPtrac);
uPtracRelative = uPtracMoy/PtracMoy*100;
EtracRef = Etrac(end);
uEtracRef = uEtrac(end);
uEtracRelative = uEtracRef/EtracRef*100;

%% Affichage

disp('g = '+string(mean(g))+' m/s² | rho = '+string(mean(rho))+' kg/m³ | Crr = '+string(mean(Crr)));
% disp(' ');
% disp('Incertitude moyenne Faero = '+string(uFaeroMoy)+' N');
% disp('Poids : rho   = '+string(prhoFaero)+' %');
% disp('        SCx   = '+string(pSCxFaero)+' %');
% disp('        v     = '+string(pvFaero)+' %');
% disp(' ');
% disp('Incertitude moyenne Frr = '+string(uFrrMoy)+' N');
% disp('Poids : Crr   = '+string(pCrrFrr)+' %');
% disp('        m     = '+string(pmFrr)+' %');
% disp('        g     = '+string(pgFrr)+' %');
% disp('        alpha = '+string(palphaFrr)+' %');
% disp(' ');
% disp('Incertitude moyenne Fgrav = '+string(uFgravMoy)+' N');
% disp('Poids : m     = '+string(pmFgrav)+' %');
% disp('        g     = '+string(pgFgrav)+' %');
% disp('        alpha = '+string(palphaFgrav)+' %');
% disp(' ');
% disp('Incertitude moyenne Fnet = '+string(uFnetMoy)+' N');
% disp('Poids : m     = '+string(pmFnet)+' % ');
% disp(' ');
disp('Proposition d incertitudes maximales pour ne pas dépasser 1% sur Etrac');
disp('        uh max = '+string(uhMax)+' m');
disp('        uT max = '+string(uTMax)+' K');
disp('        uP max = '+string(uPMax)+' P');
disp('        uSCx max = '+string(uSCxMax));
disp('        um max = '+string(umMax)+' kg');
disp('        up max = '+string(upMax)+' P');
disp('        ualpha max = '+string(ualphaMax)+' rad');
disp('        uv max = '+string(uvMax)+' m/s');
disp(' ');
disp('Part de l incertitude engendre par chaque paramètre');
disp('        h     = '+string(phP)+' %');
disp('        T     = '+string(pTP)+' %');
disp('        P     = '+string(pPP)+' %');
disp('        SCx   = '+string(pSCxP)+' %');
disp('        m     = '+string(pmP)+' %');
disp('        p     = '+string(ppP)+' %');
disp('        alpha = '+string(palphaP)+' %');
disp('        v     = '+string(pvP)+' %');
disp(' ');
disp('Ptrac moyen = '+string(PtracMoy)+' W');
disp('Incertitude moyenne Ptrac = '+string(uPtracMoy)+' W');
disp('Incertitude relative Ptrac = '+string(uPtracRelative)+' %');
disp(' ');
disp('Etrac = '+string(EtracRef/3600000)+' kWh');
disp('Incertitude Etrac = '+string(uEtracRef/3600000)+' kWh');
disp('Incertitude relative Etrac = '+string(uEtracRelative)+' %');