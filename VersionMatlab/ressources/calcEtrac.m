%Permet de calculer l'energie de traction d'un vehicule
function [Ptrac,Etrac] = calcEtrac(t,v,a,m,g,alpha,rho,S,Cx,Crr,rMot,rReg)
    %Forces
    Fgrav = m.*g.*sin(alpha);
    Faero = 1/2.*rho.*v.^2.*S.*Cx;
    Frr = Crr.*m.*g.*cos(alpha);
    Fnet = m.*a;

    Ftrac = Fgrav + Faero + Frr + Fnet;

    %Puissances
    Ptrac = Ftrac.*v;

    %Rendement en phase de motorisation et de regeneration
    for i=1:length(Ptrac)
        if(Ptrac(i)>0)
            Ptrac(i)=Ptrac(i).*(1/rMot);
        else
            Ptrac(i)=Ptrac(i).*rReg;
        end
    end
    
    %Energie en kWh
    Etrac = cumtrapz(t,Ptrac)./3600000;
end