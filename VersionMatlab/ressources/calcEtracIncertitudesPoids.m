%Permet de calculer les incertitudes des forces, de la puissance et de
%l energie, le poids des parametres qui les influence et l incertitude
%maximale de chaque parametre pour ne pas depasser un pourcentage
%d incertitude sur l energie
function [uFaero,pFaero,uFrr,pFrr,uFgrav,pFgrav,uFnet,pFnet,uPtrac,pPtrac,uEtrac,uMax,Ptrac,Etrac] = calcEtracIncertitudesPoids(t,v,a,m,alpha,S,Cx,lambda0,lambda1,lambda2,p,Crr,h,g,T,P,rhoMod,rho,uT,uP,uh,um,up,uv,ualpha,uSCx,pE,rMot,rReg)
    %% Constantes

    P0 = 101325;                    %Pression atmosphérique à h=0
    M = 0.029;                      %Masse molaire de l'air
    R = 8.314;                      %Constante des gaz parfaits
    G = 6.6742e-11;                 %Constante de gravitation universelle
    Mt = 5.9722e+24;                %Masse de la terre
    Rt = 6373000;                   %Rayon de la terre à Paris

    %% Calcul des derivees partielles et des incertitudes des parametres secondaires

    %Derivees partielles des paramètres secondaires rho, Crr et g
    drho_dP = M./(R.*T);
    drho_dT_PT = -P.*M./(R.*T.^2);
    drho_dh = -P0.*M.^2.*g./(R.^2+T.^2).*exp(-h.*M.*g./(R.*T));
    drho_dT_HT = P0.*M./(R.*T.^2).*exp(-h.*M.*g./(R.*T)).*(h.*M.*g./(R.*T)-1);

    dCrr_dp = -(lambda1 + lambda2.*v.^2)./(p.^2);
    dCrr_dv = 2.*lambda2.*v./p;

    dg_dh = (-2).*G.*Mt./(Rt+h).^3;

    %Incertitude des paramètres secondaires rho, Crr et g
    urho_PT = sqrt((drho_dT_PT.*uT).^2+(drho_dP.*uP).^2);
    urho_HT = sqrt((drho_dT_HT.*uT).^2+(drho_dh.*uh).^2);
    uCrr = sqrt((dCrr_dp.*up).^2+(dCrr_dv.*uv).^2);
    ug = sqrt((dg_dh.*uh).^2);

    %% Calcul des derivees partielles des forces

    %Derivees partielles de chaque force par rapport à ses paramètres

    %Derivees partielles de la force aerodynamique Faero
    %Faero = 1/2.*rho.*v.^2.*S.*Cx
    dFaero_drho = 0.5.*S.*Cx.*v.^2;
    dFaero_dSCx =  0.5.*rho.*v.^2;
    dFaero_dv = rho.*S.*Cx.*v;

    dFaero_dP = dFaero_drho.*drho_dP;
    dFaero_dT_PT = dFaero_drho.*drho_dT_PT;
    dFaero_dh = dFaero_drho.*drho_dh;
    dFaero_dT_HT = dFaero_drho.*drho_dT_HT;

    %Derivees partielles de la force de resistance au roulement Frr
    %Frr = Crr.*m.*g.*cos(alpha)
    dFrr_dCrr = m.*g.*cos(alpha);
    dFrr_dm = Crr.*g.*cos(alpha);
    dFrr_dg = Crr.*m.*cos(alpha);
    dFrr_dalpha = -Crr.*m.*g.*sin(alpha);

    dFrr_dp = dFrr_dCrr.*dCrr_dp;
    dFrr_dv = dFrr_dCrr.*dCrr_dv;
    dFrr_dh = dFrr_dg.*dg_dh;

    %Derivees partielles de la force gravitationnelle Fgrav
    %Fgrav = m.*g.*sin(alpha)
    dFgrav_dm = g.*sin(alpha);
    dFgrav_dg = m.*sin(alpha);
    dFgrav_dalpha = m.*g.*cos(alpha);

    dFgrav_dh = dFgrav_dg.*dg_dh;

    %Derivees partielles de la force resultante de la loi de Newton Fnet
    %Fnet = m.*a
    dFnet_dm = a;

    %% Calcul des derivees partielles de la puissance

    %Derivées partielles de la puissance par rapport à ses paramètres
    %multipliées par leurs incertitudes
    dP_drho = dFaero_drho.*v;
    dP_dSCx = dFaero_dSCx.*v;
    dP_dm = (dFrr_dm+dFgrav_dm+dFnet_dm).*v;
    dP_dg = (dFrr_dg+dFgrav_dg).*v;
    dP_dCrr = dFrr_dCrr.*v;
    dP_dalpha = (dFrr_dalpha+dFgrav_dalpha).*v;
    dP_dv = 3./2.*rho.*S.*Cx.*v.^2+(lambda0+1./p.*(lambda1+3.*lambda2.*v.^2)).*m.*g.*cos(alpha)+m.*g.*sin(alpha);

    dP_dp = dFrr_dp.*v;

    %La façon de calculer rho influence les incertitudes
    if (rhoMod == 1)
        urho = urho_HT;
        dP_dh = (dFaero_dh+dFrr_dh+dFgrav_dh).*v;
        dP_dT = dFaero_dT_HT;
        dP_dP = zeros(size(dP_dT));
    else 
        urho = urho_PT;
        dP_dh = (dFrr_dh+dFgrav_dh).*v;
        dP_dT = dFaero_dT_PT;
        dP_dP = dFaero_dP;
    end

    %% Calcul de l incertitude des forces et des poids de leurs parametres

    %Calcul de l incertitude de Faero et des poids de ses parametres
    dFaero = (dFaero_drho.*urho).^2+(dFaero_dSCx.*uSCx).^2+(dFaero_dv.*uv).^2;
    if (mean(dFaero) == 0)
        prhoFaero = 0;
        pSCxFaero = 0;
        pvFaero = 0;
    else
        prhoFaero = mean((dFaero_drho.*urho).^2)*100/mean(dFaero);
        pSCxFaero = mean((dFaero_dSCx.*uSCx).^2)*100/mean(dFaero);
        pvFaero = mean((dFaero_dv.*uv).^2)*100/mean(dFaero);
    end
    pFaero = [prhoFaero pSCxFaero pvFaero];
    uFaero = sqrt(dFaero);

    %Calcul de l incertitude de Frr et des poids de ses paramètres
    dFrr = (dFrr_dCrr.*uCrr).^2+(dFrr_dm.*um).^2+(dFrr_dg.*ug).^2+(dFrr_dalpha.*ualpha).^2;
    if (mean(dFrr) == 0)
        pCrrFrr = 0;
        pmFrr = 0;
        pgFrr = 0;
        palphaFrr = 0;
    else
        pCrrFrr = mean((dFrr_dCrr.*uCrr).^2)*100/mean(dFrr);
        pmFrr = mean((dFrr_dm.*um).^2)*100/mean(dFrr);
        pgFrr = mean((dFrr_dg.*ug).^2)*100/mean(dFrr);
        palphaFrr = mean((dFrr_dalpha.*ualpha).^2)*100/mean(dFrr);
    end
    pFrr = [pCrrFrr pmFrr pgFrr palphaFrr];
    uFrr = sqrt(dFrr);

    %Calcul de l'incertitude de Fgrav et des poids de ses paramètres
    dFgrav = (dFgrav_dm.*um).^2+(dFgrav_dg.*ug).^2+(dFgrav_dalpha.*ualpha).^2;
    if (mean(dFgrav) == 0)
        pmFgrav = 0;
        pgFgrav = 0;
        palphaFgrav = 0;
    else
        pmFgrav = mean((dFgrav_dm.*um).^2)*100/mean(dFgrav);
        pgFgrav = mean((dFgrav_dg.*ug).^2)*100/mean(dFgrav);
        palphaFgrav = mean((dFgrav_dalpha.*ualpha).^2)*100/mean(dFgrav);
    end
    pFgrav = [pmFgrav pgFgrav palphaFgrav];
    uFgrav = sqrt(dFgrav);

    %Calcul de l'incertitude de Fnet et des poids de ses paramètres
    dFnet = (dFnet_dm.*um).^2;
    if (mean(dFnet) == 0)
        pmFnet = 0;
    else
        pmFnet =100;
    end
    pFnet = [pmFnet];
    uFnet = sqrt(dFnet);

    %% Calcul de la puissance et de l energie de traction

    %Forces qui s appliquent au vehicule
    Fgrav = m.*g.*sin(alpha);       %Force gravitationnelle
    Faero = 1/2.*rho.*v.^2.*S.*Cx;  %Force aerodynamique
    Frr = Crr.*m.*g.*cos(alpha);    %Force de resistance au roulement
    Fnet = m.*a;                    %Force resultante de la loi de Newton

    %Force de traction
    Ftrac = Fgrav + Faero + Frr + Fnet;

    %Puissance de traction 
    Ptrac = Ftrac.*v;

    %Rendement en phase de motorisation et de regeneration
    for i=1:length(Ptrac)
        if(Ptrac(i)>0)
            Ptrac(i)=Ptrac(i).*(1./rMot(i));
            dP_dh(i)=dP_dh(i).*(1./rMot(i));
            dP_dT(i)=dP_dT(i).*(1./rMot(i));
            dP_dP(i)=dP_dP(i).*(1./rMot(i));
            dP_dp(i)=dP_dp(i).*(1./rMot(i));
            dP_dSCx(i)=dP_dSCx(i).*(1./rMot(i));
            dP_dm(i)=dP_dm(i).*(1./rMot(i));
            dP_dalpha(i)=dP_dalpha(i).*(1./rMot(i));
            dP_dv(i)=dP_dv(i).*(1./rMot(i));
        else
            Ptrac(i)=Ptrac(i).*rReg(i);
            dP_dh(i)=dP_dh(i).*rReg(i);
            dP_dT(i)=dP_dT(i).*rReg(i);
            dP_dP(i)=dP_dP(i).*rReg(i);
            dP_dp(i)=dP_dp(i).*rReg(i);
            dP_dSCx(i)=dP_dSCx(i).*rReg(i);
            dP_dm(i)=dP_dm(i).*rReg(i);
            dP_dalpha(i)=dP_dalpha(i).*rReg(i);
            dP_dv(i)=dP_dv(i).*rReg(i);
        end
    end

    %Energie en kWh
    Etrac = cumtrapz(t,Ptrac);

    %% Calcul de l incertitude de la puissance et des poids de ses parametres

    %Calcul de l'incertitude de Ptrac et des poids de ses paramètres
    dPtrac = (dP_dh.*uh).^2+(dP_dT.*uT).^2+(dP_dP.*uP).^2+(dP_dp.*up).^2+(dP_dSCx.*uSCx).^2+(dP_dm.*um).^2+(dP_dalpha.*ualpha).^2+(dP_dv.*uv).^2;
    if (any(dPtrac))
        phP = mean((dP_dh.*uh).^2)*100/mean(dPtrac);
        pTP = mean((dP_dT.*uT).^2)*100/mean(dPtrac);
        pPP = mean((dP_dP.*uP).^2)*100/mean(dPtrac);
        pSCxP = mean((dP_dSCx.*uSCx).^2)*100/mean(dPtrac);
        pmP = mean((dP_dm.*um).^2)*100/mean(dPtrac);
        ppP = mean((dP_dp.*up).^2)*100/mean(dPtrac);
        palphaP = mean((dP_dalpha.*ualpha).^2)*100/mean(dPtrac);
        pvP = mean((dP_dv.*uv).^2)*100/mean(dPtrac);
    else
        phP = 0;
        pTP = 0;
        pPP = 0;
        pSCxP = 0;
        pmP = 0;
        ppP = 0;
        palphaP = 0;
        pvP = 0;
    end
    pPtrac = [phP pTP pPP pSCxP pmP ppP palphaP pvP];
    uPtrac = sqrt(dPtrac);

    %Calcul de l'incertitude de Etrac
    uEtrac = cumtrapz(t,uPtrac);

    %% Calcul de l incertitude de chaque parametre pour ne pas depasser pE% sur Etrac
    
    %Nombre de parametres pris en compte
    if (rhoMod == 1)
        nP = sum([uh uT uSCx um up ualpha uv]==[0 0 0 0 0 0 0]);
    else
        nP = sum([uh uT uP uSCx um up ualpha uv]==[0 0 0 0 0 0 0 0]);
    end
    
    %Energie de traction restante a repartir sur les autre parametre
    poids = pE/100*Etrac(end)/sqrt(nP);
    if(any(dP_dh)&&(uh==0))
        ph = cumtrapz(t,abs(dP_dh));
        uhMax = poids./ph(end);
    else
        uhMax = 0;
    end
    if(any(dP_dT)&&(uT==0))
        pT = cumtrapz(t,abs(dP_dT));
        uTMax = poids./pT(end);
    else
        uTMax = 0;
    end
    if(any(dP_dP)&&(uP==0))
        pP = cumtrapz(t,abs(dP_dP));
        uPMax = poids./pP(end);
    else
        uPMax = 0;
    end
    if(any(dP_dSCx)&&(uSCx==0))
        pSCx = cumtrapz(t,abs(dP_dSCx));
        uSCxMax = poids./pSCx(end);
    else
        uSCxMax = 0;
    end
    if(any(dP_dm)&&(um==0))
        pm = cumtrapz(t,abs(dP_dm));
        umMax = poids./pm(end);
    else
        umMax = 0;
    end
    if(any(dP_dp)&&(up==0))
        pp = cumtrapz(t,abs(dP_dp));
        upMax = poids./pp(end);
    else
        upMax = 0;
    end
    if(any(dP_dalpha)&&(ualpha==0))
        palpha = cumtrapz(t,abs(dP_dalpha));
        ualphaMax = poids./palpha(end);
    else
        ualphaMax = 0;
    end
    if(any(dP_dv)&&(uv==0))
        pv = cumtrapz(t,abs(dP_dv));
        uvMax = poids./pv(end);
    else
        uvMax = 0;
    end

    uMax = [uhMax uTMax uPMax uSCxMax umMax upMax ualphaMax uvMax];
end