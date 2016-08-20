function dydt=ODE_LHS(t,y,LHSmatrix,x,runs)
%% PARAMETERS %%
Parameter_settings_LHS;
mu=LHSmatrix(x,1);
gamma=LHSmatrix(x,2);
sigma=LHSmatrix(x,3);
omega=LHSmatrix(x,4);
rho=LHSmatrix(x,5);
psi=LHSmatrix(x,6);
beta=LHSmatrix(x,7);
% #s=LHSmatrix(x,1);
% #muT=LHSmatrix(x,2);
% #r=LHSmatrix(x,3);
% #k1=LHSmatrix(x,4);
% #k2=LHSmatrix(x,5);
% #mub=LHSmatrix(x,6);
% #N=LHSmatrix(x,7);
% #muV=LHSmatrix(x,8);
dummy_LHS=LHSmatrix(x,8);

    dydt(1) = - beta * y(1) * y(2)./(y(1)+y(2)+y(3)+y(4))- mu * y(1); 
    dydt(2) = beta*y(1) * y(2)./(y(1)+y(2)+y(3)+y(4)) +gamma * y(3) - mu * y(2)- sigma*y(2)-psi*y(2) + omega*y(4); 
    dydt(3) = sigma*y(2) -gamma * y(3) - mu * y(3)-rho*y(3); 
    dydt(4) = rho*y(3) - mu*y(4) -omega*y(4); 
    dydt=dydt';
    
%     ##############################################################################################################3

% #% [T] CD4+ uninfected: Tsource + Tprolif - Tinf
% #Tsource = s - muT*y(1);
% #Tprolif = r*y(1)*(1-(y(1)+y(2)+y(3))/Tmax);
% #Tinf = k1*y(1)*y(4);
% 
% #% [T1] CD4+ latently infected: Tinf - T1death - T1inf
% #T1death = muT*y(2);
% #T1inf = k2*y(2);
% 
% #% [T2] CD4+ actively infected: T1inf - T2death
% #T2death = mub*y(3);
% 
% #% [V] Free infectious virus: Vrelease - Tinf - Vdeath
% #Vrelease = N*T2death;
% #Vdeath = muV*y(4);
% 
% #dydt = [Tsource + Tprolif - Tinf;
% #        Tinf - T1death - T1inf;
% #        T1inf - T2death;
% #        Vrelease - Tinf - Vdeath];
