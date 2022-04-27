%Retourne le temps, la distance, la vitesse et l'acceleration du cycle WLTP
function [temps,distance,vitesse,acceleration] = getWLTP()

WLTP = load('WLTP.mat');

temps = WLTP.WLTP_class3(:,1);
vitesse = WLTP.WLTP_class3(:,2)/3.6;
tempsv=diff(temps);
tempsv=[1; tempsv];
acceleration=diff(vitesse);
acceleration=[0;acceleration];
acceleration=acceleration./tempsv;
distance = cumtrapz(temps,vitesse);
end