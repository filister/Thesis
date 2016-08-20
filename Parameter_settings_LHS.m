% PARAMETER BASELINE VALUES

%  mu=0.0615./2;
%  gamma=0.0486./2;
% sigma=0.0027./2;
%  omega=3.6807e-007./2;
%  rho=0.0379./2;
% psi=0.0042./2;
% beta=0.1547./2;
%age 12
% mu=6.5587e-007;
% gamma=0.0170;
% sigma=0.0393;
% omega=0.0298;
% rho=1.2571e-006;
% psi=4.5216e-006;
% beta=0.2694;

%age 22
mu=9.5647e-009 ;
gamma=0.0290;
sigma=0.1554;
omega=0.0853;
rho=7.6910e-009;
psi=0.0031 ;
beta=6.3333e-007;
%age 17
% mu=9.5135e-009;
% gamma=8.2820e-006;
% sigma=0.2680;
% omega=0.1196;
% rho=0.2894;
% psi=9.5789e-007;
% beta=0.1713;

% mu=0.0037./1000;
% gamma=1./1000;
% sigma=0.0059;
% omega=0.1058;
% rho=0.8599;
% psi=0.0085;
% beta=0.0134;
% #s=10; 
% #muT=2e-2;
% #r=0.03;
% #Tmax=1500;
% #k1=2.4e-5;
% #k2=3e-3;
% #mub=0.24;
% #N=1200;
% #muV=2.4;
%age 27
% mu=3.0221e-007;
% gamma=0.2276;
% sigma=0.0017;
% omega=0.3613;
% rho=0.2196;
% psi=0.0532;
% beta=0.1819;

mu=3.0221e-007;
gamma=0.2276;
sigma=0.0017;
omega=0.3613;
rho=0.2196;
psi=0.0532;
beta=0.1819;
dummy=1;

% Parameter Labels 
PRCC_var={'mu', 'gamma', 'sigma', ...
    'omega','rho', 'psi','beta','dummy'};% 

%% TIME SPAN OF THE SIMULATION
t_end=100; %length of the simulations
tspan=(0:1:t_end);   % time points where the output is calculated
time_points=[10,t_end./2,t_end]; % time points of interest for the US analysis

% INITIAL CONDITION FOR THE ODE MODEL
S0=68568; D0=18589; R0=115.2; Q0=10;
% #T0=1e3;
% #T1=0;
% #T2=0;
% #V=1e-3;
% #y0=[T0,T1,T2,V];
y0=[S0,D0,R0,Q0];

% Variables Labels
y_var_label={'S','D','R','Q'};
