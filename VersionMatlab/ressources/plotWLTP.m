%% Plot WLTP class 3b

clear;
close all;
clc;

addpath(genpath('../etude_sensibilite'));

[t,d,v,a] = getWLTP();


figure(1)
plot(t,v)
title("WLTP class 3b - Profil de vitesse")
xlabel('Temps (s)')
ylabel('Vitesse (m/s)')
xline(589)
xline(1022)
xline(1477)
