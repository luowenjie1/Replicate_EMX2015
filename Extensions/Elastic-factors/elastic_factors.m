clear;
clc;


% neoclassical growth model with elastic labor supply

alpha = 1/3 ; %% output elasticity of capital in aggregate production function
delta = 0.1;  %% capital depreciation
beta  = 0.96; %% discount factor
L     = 1;    %% benchmark labor supply, choose chi to ensure this   
eta   = 1;    %% inverse Frisch elasticity of labor supply

rho   = (1/beta) - 1; %% pure rate of time preference


% aggrgate statistics from benchmark model

A_autarky     = 23.204;    %% aggregate productivity in autarky equilibrium
A_taiwan      = 26.274;    %% aggregate productivity in benchmark taiwan calibration

A_autarky_eff = 25.377;    %% efficient productivity in autarky equilibrium
A_taiwan_eff  = 28.1722;   %% efficient productivity in benchmark taiwan calibration

mu_autarky    = 1.353;     %% aggregate markup       in autarky equilibrium
mu_taiwan     = 1.314;     %% aggrgeate markup       in benchmark taiwan calibration


% notes:

% -- 10.4% increase in productivity from A_autarky_eff to A_taiwan_eff
% -- 12.4% increase in productivity from A_autarky     to A_taiwan
% --  2.9% decrease in markup       from mu_taiwan     to mu_autarky


%%%%% compute terminal steady state

Abar = A_taiwan;
Mbar = mu_taiwan;

chi  = (1-alpha)/(Mbar - (alpha*delta)/(rho+delta)); %% weight on leisure that ensures Lbar = 1

parameters = [alpha,delta,beta,chi,eta];

[Cbar,Kbar,Ybar,Lbar,~] = elastic_factors_steady_state(parameters,Abar,Mbar);


%%%%% compute initial steady state

A0 = A_autarky;
M0 = mu_autarky;


[C0,K0,Y0,L0,~] = elastic_factors_steady_state(parameters,A0,M0);

%%%% compute transitional dynamics from initial to terminal steady state

T  = 500; 

Cguess = ((2/3)*Cbar+(1/3)*C0)*ones(T,1);
Kguess = ((2/3)*Kbar+(1/3)*K0)*ones(T,1); %% initial guess

Xguess = [Cguess,Kguess];
XLB    = zeros(T,2);
XUB    = inf(T,2);

options = optimoptions('lsqnonlin');
options = optimoptions(options,'MaxFunEvals',1000*T);

Xsolution = lsqnonlin(@(X) elastic_factors_transitional_dynamics(X,parameters,Abar,Mbar,Kbar,K0,Cbar,Lbar),Xguess,XLB,XUB,options);


Cpath = Xsolution(:,1);
Kpath = Xsolution(:,2);

Lpath = (((1-alpha)/chi)*(Abar/Mbar)*(Kpath.^(alpha))./Cpath).^(1/(eta+alpha)); %% closed-form solution for L given C,K

Ypath = Abar*(Kpath.^alpha).*(Lpath.^(1-alpha));

Upath = log(Cpath)-chi*(Lpath.^(1+eta))/(1+eta);

U1    = sum((beta.^(0:1:T-1))'.*Upath);            %% utility from time t=0, including transition path
U0    = (log(C0)-chi*L0^(1+eta)/(1+eta))/(1-beta); %% utility associated with staying at initial steady state~

fprintf('\n')
display('long-run % change in C,K,Y,L from given % change in A,M')
fprintf('\n')
display('capital accumulation, elastic labor supply')
fprintf('\n')
fprintf('A change, percent                      = %7.3f \n', log(Abar/A0)*100);
fprintf('M change, percent                      = %7.3f \n', log(Mbar/M0)*100);
fprintf('\n')
fprintf('C change, percent                      = %7.3f \n', log(Cbar/C0)*100);
fprintf('K change, percent                      = %7.3f \n', log(Kbar/K0)*100);
fprintf('Y change, percent                      = %7.3f \n', log(Ybar/Y0)*100);
fprintf('L change, percent                      = %7.3f \n', log(Lbar/L0)*100);
fprintf('\n')
fprintf('compensating variation, percent        = %7.3f \n', log(exp((1-beta)*(U1-U0)))*100);
fprintf('(including transition)')
fprintf('\n')
fprintf('\n')



