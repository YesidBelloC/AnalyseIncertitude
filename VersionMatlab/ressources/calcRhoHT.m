%Permet de calculer la masse volumique de l'air en fonction de l'altitude
function rho = calcRhoHT(h,T,g)
    P0 = 101325;
    M = 0.029;
    R = 8.314;
    rho = P0.*M.*exp(-h.*M.*g./(R.*T))./(R.*T);
end