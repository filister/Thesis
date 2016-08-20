

function drugabusefitw(do_estimation)
warning off;

P_Data(:,1)=[1:19];
P_Data(:,2)=[1.6 3.8 3.6 2.6 4.2 3.2 2.8 2.8 1.9 2.2 2.2 2.2 1.6 2.6 2.4 3.8 4.2 2.2 3.8]./100 ;%age 67 data
Tspan=P_Data(:,1);
% Parameters
beta= 0.3 ;
mu =5e-2;
gamma =0.2;
sigma =0.4;
omega =0.3;
rho =0.3;
psi =0.5;
   
S0=20919; D0=57; R0=1.6; Q0=0.3;  
INITIAL=[S0,D0,R0,Q0]/100;
Istart=1;    % month to start the model simulation
Iend=Istart+ 19;

OPTIONS=odeset('AbsTol',0.001,'RelTol',0.001,'MaxStep',1/12);
% tdur=50;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate parameters 
% by minimizing the sum of squares
% when fitting modeled to real prevalence data
    
do_estimation=1;
if(do_estimation)
xdata=P_Data(:,1)';
ydata=P_Data(:,2)';
x0(1,1)= 0.13;       % beta;
x0(1,2)=5e-2;      % mu;
x0(1,3)= 0.02;  % gamma;
x0(1,4)= 0.0045 ;     % sigma;
x0(1,5)=0.01  ;    % omega;
x0(1,6)=0.0067 ;       % rho;
x0(1,7)=0.000001  ;    %psi;
% LB=[0 0.0 0.1 0.1 0.0 0.0 0.0 0.0 0.0 10 0.0 0.0 0.0 0.0 0.0 0.0];
% UB=[1 0.2 0.9 0.9 0.8 1 0.18 0.5 0.7 20 1 0.1 0.3 0.7 0.9 0.1];
LB=[0 0 0 0 0 0 0 ];
UB=[ 1 1 1 1 1 1 1 ];
 options = optimset('MaxFunEvals', 100000,'MaxIter',45000,'TolFun',1e-9,'TolX',1e-9);
% LB(12)=UB(12)
% options = optimset('MaxFunEvals',10000);
%options = odeset('RelTol',1e-4,'AbsTol',1e-50*ones(1,4));
x=lsqcurvefit(@Model_Inc,x0,xdata,ydata,LB,UB,options);
% x=lsqnonlin(@Model_Inc,x0,xdata,ydata,LB,UB,options);
% x = fit(@Model_Inc,xdata,ydata,x0,LB,UB,'Robust','on');

'estimated parameters'
beta= x(1)
mu =x(2)
gamma =x(3)
sigma =x(4)
omega =x(5)
rho =x(6)
psi =x(7)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[t y] = ode45(@drugabusefitw,[0:1:(Iend-Istart)], INITIAL);
S=y(:,1);
D=y(:,2);
R=y(:,3);
Q=y(:,4);

%incidence= (y(:,2)*y(:,1)*beta+ y(:,3)*gamma +  y(:,4)*omega )/(y(:,1) + y(:,3) + y(:,4) );
incidence= y(:,3);
% incidence= y(:,4)*gamma;
close all;
figure(1);
hold on
h_l=plot(Istart+t,incidence,'r-');
set(h_l,'linewidth',2);
h_l=plot(P_Data(:,1),P_Data(:,2),'bo','Markersize',8);
set(h_l,'linewidth',2);
axis([Istart Iend 0 0.45]);
ylim([0 0.45]) %                              The line added to change the y limit
% title('Rehabilitation Prevalence (%)');
xlabel('years','fontsize',15)
ylabel('Drug abuse Incidence','fontsize', 15)
hold off

function [dydt]=drugabusefitw(t,y)

S=y(1);
D=y(2);
R=y(3);
Q=y(4);


    dydt(1) = - beta * y(1) * y(2)./(y(1)+y(2)+y(3)+y(4))- mu * y(1); 
    dydt(2) = beta*y(1) * y(2)./(y(1)+y(2)+y(3)+y(4)) +gamma * y(3) - mu * y(2)- sigma*y(2)-psi*y(2) + omega*y(4); 
    dydt(3) = sigma*y(2) -gamma * y(3) - mu * y(3)-rho*y(3); 
    dydt(4) = rho*y(3) - mu*y(4) -omega*y(4); 
    dydt=dydt';
 
end % function 


function Inc=Model_Inc(x0,xdata)

Inc=0;      % intialization of this not to have an empty array
beta=x0(1);
mu=x0(2);
gamma=x0(3);
sigma=x0(4);
omega=x0(5);
rho=x0(6);
psi=x0(7);
%%% Initial values 
% S0=35000; E0=0; Ia0=0; Is0=1; R0=0;
S0=20919; D0=57; R0=1.6; Q0=0.3; 
INITIAL=[S0,D0,R0,Q0]/100;
OPTIONS=odeset('AbsTol',0.001,'RelTol',0.001,'MaxStep',1/2);
tdur=50;
[t y] = ode45(@drugabusefitw, [0:1:(Iend-Istart)], INITIAL);

S=y(:,1);
D=y(:,2);
R=y(:,3);
Q=y(:,4);

  

%incidence= (1-p)*gamma*y(:,2)+ theta*y(:,3)
% incidence= y(:,4)*gamma;

for i=1:length(xdata)
    
    ind=find(xdata(i)-Istart==t);
     S=y(ind,1);
     D=y(ind,2);
     R=y(ind,3);
     Q=y(ind,4);
            Inc(i)= R;  % incidence
     % Inc(i)= Is*gamma;  % incidence
       % Inc(i)= (beta*D*S + gamma*R +  omega*Q )./(S + R + Q )
end 

end     % end of Model_prev    

% Q1= gamma + mu; 
% Q2= mu + theta + delta1;
% Q3 = delta2 + sigma + mu;
% mstar= Lambda*alpha1 / omega* mu;
% R_M = beta*c*gamma*(1-mstar)*(p*theta + Q2*(1-p))/(Q1*Q2*Q3)
 end     % end of ebola_fit
  
