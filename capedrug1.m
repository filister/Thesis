clear all

h = 0.2;  % h is the grid size
begin_year = 2005.5;
age = 12:h:67;
time = begin_year:h:2020;
[X,Y] = meshgrid(age,time);

% Variable which runs the code for a number of times to provide convergence
numruns = 100;

% Declare all the classes and incidence by means of zero matrices 
S = zeros(length(age),length(time));
D = zeros(length(age),length(time));
R = zeros(length(age),length(time));
Q = zeros(length(age),length(time));
Incidence = zeros(length(age),length(time));


Proportion = 23164628/1312480000;
% Boundary Condition
S(1,:) = 90180*exp(0.0261*(time-2005.5)); 

% Initial Conditions
%S(:,1) = [700 500 580 690 650 543 488.2 401 350 164.8 259.6 200];#column 1 of the matrix of susceptibles
%D(:,1) = [100 200 190 75 90 45.2 63.6 30 25 50 20 13];  
%R(:,1) = [19.6 115.2 101 51.4 40 31.8 27.4 17.8 12.2 5.4 2 1.6]; 
%Q(:,1) = [7 15 20 30 40 20 5 50 7.2 100 8 0.4]; 

% S(:,1)=635.26 +(2.071-0.15*age).*(age);
% D(:,1)=222.19 +(0.03*age-5.12).*(age);
% R(:,1)=87-(1.01+0.006*age).*(age);
% Q(:,1)=2.61*age-(0.03*age.^2 + 21.8);

S(:,1)=61616 +(258-13*age).*(age);
%D(:,1)=24328 -(6*age-21).*(age);
D(:,1)=222.19 +(0.03*age-5.12).*(age);
R(:,1)=87-(1.01+0.006*age).*(age);
Q(:,1)=17.4028-(0.0013*age + 0.2022).*(age);
%Code to find right population size and proportion of classes 
sum(S(:,1) + D(:,1)+ D(:,1)+ Q(:,1))
[sum(S(:,1)) sum(D(:,1)) sum(R(:,1)) sum(Q(:,1)) ]
   sum(S(:,1) + D(:,1)+ R(:,1)+ Q(:,1))


% Ages linked to the corresponding entries of parameter vectors 
ages = [12 17 22 27 32 37 42 47 52 57 62 67];

% Parameters vectors corresponding to ages
mu=[0.0615 0.0037 3.5735e-009 3.0221e-007 1.4811e-007 1.1365e-008 6.0586e-008 0.0011 0.0091 0.0048 3.9316e-005 3.8089e-008];
gamma=[0.0486 1 0.4571 0.2276 0.9878 1 0.1021 0.0798 0.2905 0.7962 1 0.99];
sigma=[0.0027 0.0059 0.0031 0.0017 0.0063 0.0080 0.0043 0.0034 0.0036 0.0098 0.0156 0.1271];
omega=[3.6807e-007 0.1058 0.4344 0.3613 0.2182 0.1562 0.1055 0.1254 0.1119 0.1178 0.0898 0.1209];
rho=[0.0379 0.8599 0.4404 0.2196 0.99 1 0.9644 0.9898 0.6315 0.9676 1 0.9851];
psi=[0.0042 0.0085 0.0714 0.0532 2.3094e-014 7.0021e-004 0.0089 0.0147 0.0069 0.0229 0.222 0.1075];
% mu=[6.5587e-007 9.5135e-009 9.5647e-009 1.1919e-005 1.0193e-009 4.3636e-005 4.3429e-004 0.0524 5.4327e-004 0.0286 8.0360e-004 2.2623e-004];
% gamma=[0.0170 8.2820e-006 0.0290 6.3804e-013 2.8180e-008 0.0521 0.0164 5.2131e-004 0.0540 0.2006 0.0047 7.1540e-005];
% sigma=[0.0393 0.2680 0.1554 0.0955 0.0345 0.7873 0.9357 0.7354 0.3090 0.3479 0.9722 7.5140e-008];
% omega=[0.0248 0.1196 0.0853 0.0326 0.0253 0.0045 0.0873 0.9931 8.1030e-004 2.9466e-004 0.0706 0.7431];
% rho=[1.2571e-006 0.2894 7.6910e-009 2.3781e-009 0.0128 0.0015 0.1806 0.3478 2.4186e-007 9.0851e-007 0.6019 0.6445];
% psi=[4.5216e-006 9.5789e-007 0.0031 0.0097 0.0069 3.6722e-007 1.5630e-007 1.8977e-007 1.8195e-007 9.0948e-013 7.0783e-012 0.0036];


% Declare all the transmission matrix by means of a zero matrix 
betamatty = zeros(length(time),length(age));
years = [ 2005:1:2014 2020];  

% Specify transmission vectors which entries correspond to:

% Ages      12 17 22 27 32 37 42 47 52 57 62 67

beta_2004 = [0.01547 0.134 0.1726 0.1819 0.0642 0.0163 0.0128 0.0126 0.0258 0.0330 0.02164 0.01176];
beta_2005 = [0.01547 0.134 0.1726 0.1819 0.0642 0.0163 0.0128 0.0126 0.0258 0.0330 0.02164 0.01176];
beta_2006 = [0.01547 0.134 0.1726 0.1819 0.0642 0.0163 0.0128 0.0126 0.0258 0.0330 0.02164 0.01176];
beta_2007 = [0.01547 0.134 0.1726 0.1819 0.0642 0.0163 0.0128 0.0126 0.0258 0.0330 0.02164 0.01176];
beta_2008 = [0.01547 0.134 0.1726 0.1819 0.0642 0.0163 0.0128 0.0126 0.0258 0.0330 0.02164 0.01176];
beta_2009 = [0.01547 0.134 0.1726 0.1819 0.0642 0.0163 0.0128 0.0126 0.0258 0.0330 0.02164 0.01176];
beta_2010 = [0.01547 0.134 0.1726 0.1819 0.0642 0.0163 0.0128 0.0126 0.0258 0.0330 0.02164 0.01176];
beta_2011 = [0.01547 0.134 0.1726 0.1819 0.0642 0.0163 0.0128 0.0126 0.0258 0.0330 0.02164 0.01176];
beta_2012 = [0.01547 0.134 0.1726 0.1819 0.0642 0.0163 0.0128 0.0126 0.0258 0.0330 0.02164 0.01176];
beta_2013 = [0.01547 0.134 0.1726 0.1819 0.0642 0.0163 0.0128 0.0126 0.0258 0.0330 0.02164 0.01176];
beta_2014 = [0.01547 0.134 0.1726 0.1819 0.0642 0.0163 0.0128 0.0126 0.0258 0.0330 0.02164 0.01176];
beta_2020 = [0.01547 0.134 0.1726 0.1819 0.0642 0.0163 0.0128 0.0126 0.0258 0.0330 0.02164 0.01176];
% Assigning the parameter vectors to their respective matrices

bmatty = vertcat(beta_2005,beta_2006,beta_2007,beta_2008,beta_2009,beta_2010,beta_2011,beta_2012,beta_2013,beta_2014,beta_2020)/10000; 

paramatty = vertcat(mu, gamma, sigma, omega, rho, psi);



% Code for interpolating parameters linearly between different ages 
size_paramatty = size(paramatty);
nrows_paramatty = size_paramatty(1);
paramdiffy = zeros(nrows_paramatty,length(ages)-1);
param = zeros(nrows_paramatty,length(age));

for i = 1:1:nrows_paramatty
    for j = 1:1:(length(ages)-1)
        paramdiffy(i,j) = ((paramatty(i,j+1) - paramatty(i,j))*h)/(ages(j+1) - ages(j));
    end
end    

for k = 1:1:(nrows_paramatty)
    vect1 = [];
    for m = 1:1:(length(ages)-2)
        vect1 = horzcat(vect1,paramatty(k,m):paramdiffy(k,m):(paramatty(k,m+1)-paramdiffy(k,m)));
    end
    vect1 = horzcat(vect1,paramatty(k,(length(ages)-1)):paramdiffy(k,(length(ages)-1)):paramatty(k,length(ages)));
    param(k,:) = vect1; 
end
mu = param(1,:);
gamma = param(2,:);
sigma= param(3,:);
omega = param(4,:);
rho = param(5,:);
psi = param(6,:);


ages = [12 17 22 27 32 37 42 47 52 57 62 67];
bdiffy = zeros(length(years),length(ages)-1);

for i = 1:1:length(years)
    for j = 1:1:(length(ages)-1)
        bdiffy(i,j) = ((bmatty(i,j+1) - bmatty(i,j))*h)/(ages(j+1) - ages(j));
    end
end    
  
% Variable for converting a year number into a year since the initial year
yearsad = years - begin_year; 

% Code for interpolating transmission parameters linearly between different ages and times
for k = 1:1:length(years)
    vect2 = [];
    for m = 1:1:(length(ages)-2)
        vect2 = horzcat(vect2,bmatty(k,m):bdiffy(k,m):(bmatty(k,m+1)-bdiffy(k,m)));      
    end
    vect2 = horzcat(vect2,bmatty(k,(length(ages)-1)):bdiffy(k,(length(ages)-1)):bmatty(k,length(ages)));
    betamatty=kron(ones(81,1),vect2)'; 
end
%Double for loop of the Infectives Class'es first guess.
for j = 2:1:length(time)       
    for i = 2:1:length(age)  
        lambda = betamatty(i-1,j)*D(i-1,j);
        S(i,j) = (S(i-1,j)+S(i,j-1))/(2+h*(lambda + mu(i-1)));
        D(i,j)=(D(i-1,j)+D(i,j-1)+ h*(lambda*S(i,j)+ gamma(i-1)*R(i,j)+ omega(i-1)*Q(i,j)))/ (2+ h*(mu(i-1) +sigma(i-1) + psi(i-1)));
        R(i,j) = (R(i-1,j)+R(i,j-1)+h*sigma(i-1)*D(i,j))/(2+h*(mu(i)+rho(i-1)+gamma(i-1)));
        Q(i,j)=(Q(i-1,j)+Q(i,j-1)+h*rho(i-1)*R(i,j))/(2+h*(mu(i-1) +omega(i-1)));
    end       
end

% Triple for loop of the convergence for the system of difference equations.
for runs = 1:1:numruns
    for j = 2:1:length(time) 
        for i = 2:1:length(age)
        lambda = betamatty(i,j)*D(i,j);
        S(i,j) = (S(i-1,j)+S(i,j-1))/(2+h*(lambda + mu(i)));
        D(i,j)=(D(i-1,j)+D(i,j-1)+ h*(lambda*S(i,j)+ gamma(i)*R(i,j)+ omega(i)*Q(i,j)))/ (2+ h*(mu(i) +sigma(i) + psi(i)));
        R(i,j) = (R(i-1,j)+R(i,j-1)+h*sigma(i)*D(i,j))/(2+h*(mu(i)+rho(i)+gamma(i)));
        Q(i,j)=(Q(i-1,j)+Q(i,j-1)+h*rho(i)*R(i,j))/(2+h*(mu(i) +omega(i)));
        end
    end
end

% Double for loop to fill out the matrix of incidence.
for j = 2:1:length(time)
    for i = 2:1:length(age)            
        Incidence(i,j) = (betamatty(i,j)*D(i,j)*S(i,j)); 

    end        
end      

% Years and ages related to the Cape Town Metropole Data
years = [2005.5:0.5:2014.5];

ages = [12 17 22 27 32 37 42 47 52 57 62 67];
U = zeros(length(years),length(ages));

% Cape Town Metropole Data
U(1,1:12) = [ 98  576 505 257 200 159  137  89  61 27  10 8];
U(2,1:12) = [77 647 579 340 284 239 174  134 80 51 27  19]; 
U(3,1:12) = [108 652 683 342 279 241 182 143 75 47 25 18];
U(4,1:12) = [112 688 661 357 253 230 226 149 85 56  25 13];
U(5,1:12) = [107 704 708 437 276 254 201 157 84 76 25 21];
U(6,1:12) = [71 553 582 435 248 214 172 150 97 55 36 16];
U(7,1:12) = [85 569 691 463 278 193 187 149 91 53 23 14];
U(8,1:12) = [128 764 848 633 358 324 237 157 103 62 25 14];
U(9,1:12) = [61 545 607 488 283 215 179 120 60 42 21 9]; 
U(10,1:12) = [98 595 676 627 339 258 216 114 92 41 28 11];
U(11,1:12) = [112  688 661 357 253 230 226 149 85  56  25  13];
U(12,1:12) = [119 691 676 576 334 245 193 144  92 41 28 11];
U(13,1:12) = [76 353 694 605 334 205 166 137 79 44 16 8];
U(14,1:12) = [170 696 886 845 466 276 196 151 100 56 25 13];
U(15,1:12) = [123 531 629 674 433 255 166 155 90 53 28 12]; 
U(16,1:12) = [125 617 751 825 489 308 223 146 104 62 24 1];
U(17,1:12) = [187 701 574 755 459 264 195 140 85 62 19 21];
U(18,1:12) = [167 635 561 796 504 255 219 159 96 65 18 11];
U(19,1:12) = [185 597 561 725 501 302 201 154 92 53 28 19];

% Length of age intervals
int_length = [5*ones(1,12)];

% Divide each data entry by its corresponding age interval
for i = 1:1:length(years)
      U(i,1:12) =  U(i,1:12)./int_length;     
end


% Code for scatterplot
yearsmat = zeros(length(years),length(ages));
for i = 1:1:length(years)
   yearsmat(i,:) = ones(1,length(ages))*years(i);
end

Uvec = [];
agesvec = [];
yearsvec = [];
for i = 1:1:length(years)
    Uvec = horzcat(Uvec,U(i,:));
    agesvec = horzcat(agesvec,ages);
    yearsvec = horzcat(yearsvec,yearsmat(i,:));    
end

% Code for only a fit to the data and no projections are considered. 
Xadj1 = X((5.4/h+1):1:length(time),:);
Yadj1 = Y(1:1:(9/h+1),:);
IncidenceAdj1 = Incidence(:,1:1:(9/h+1));

% Code for a fit to the data together with projections.
Xadj2 = X(1:1:length(time),:);
Yadj2 = Y(1:1:length(time),:);
IncidenceAdj2 =D(:,1:1:length(time));

figure(1)
hold on
%mesh(Xadj1,Yadj1,IncidenceAdj1')  % Alternative to only considering data
mesh(Xadj2,Yadj2,IncidenceAdj2') % Plot for the projection
scatter3(agesvec,yearsvec,Uvec,25,'k','filled')   
xlabel('age')
ylabel('time')
zlabel('Incidence')
hold off
%save capedrug1.mat
   
