clear all;
clc;

global gamma theta tau ebar N S

% Assigned Parameters

gamma   = 10.4971;                    
theta   = 1.2351;

S       = 1000;                 % number of industries: must be even
L       = 1;                    % inelastic labor supply in Home
Ls      = 1;                    % inelastic labor supply in Foreign

lambda  = 0.15;                 % prob. all large guys collude

% Calibrated

N       = 140;                                  % number of firms that draw
ebar    = 1.7477;                          % productivity of type 2
su      = 0.1706;                         % standard deviation of Home/Foreign productivity gap
tau     = 1.1269; 

om1     = [0.1085 ; 0.0128;  0.3006]; 
om2     = [0.0642; 0.0976; 0.0516]; 

pom     = [0.9147; 0.5660;  0.2033; ...
           0.4803; 0.4149;  0.4887; ...
           0.5512; 0.6701;  0.4749];
  
pom      = pom/sum(pom); 
       
[uu, wu] = qnwunif(7, 0, 1);
uu       = exp(norminv(uu, 0, su)); 

equilibrium

% % Simple example
% 
% gamma = 0.5; 
% N     = 2; 
% omega = 1/3; 
% n     = (1:1:N)';
% 
% pn    = binopdf(n, N, omega);
% 
% Aprofits = pn'*(n.*n.^(-gamma))/N; 
% Eprofits = omega*pn'*n.^(-gamma); 