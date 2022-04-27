%Permet de calculer le coefficient de résistance au roulement
function Crr = calcCrr(lambda0,lambda1,lambda2,p,v)
    Crr =lambda0+1./p.*(lambda1+lambda2.*v.^2);
end