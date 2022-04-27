%Permet de calculer l'intensité gravitationnelle
function g = calcG(h)
    G = 6.6742e-11;                 %Constante de gravitation universelle
    Mt = 5.9722e+24;                %Masse de la terre
    Rt = 6373000;                   %Rayon de la terre à Paris
    g = (G*Mt)./((Rt+h).^2);          
end