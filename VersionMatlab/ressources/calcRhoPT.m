%Permet de calculer la masse volumique de l'air en fonction de la pression
%atmosphérique
function rho = calcRhoPT(P,T)
    M = 0.029;
    R = 8.314;
    rho = P.*M./(R.*T);
end