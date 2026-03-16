clear all;
clc;

global gamma theta tau ebar N S

% Assigned Parameters

gamma   = 10.4971;                    
theta   = 1.2351;

S       = 1000;                 % number of industries: must be even
L       = 1;                    % inelastic labor supply in Home
Ls      = 1;                    % inelastic labor supply in Foreign

lambda  = 0.35;                 % prob. all large guys collude

% Calibrated

N       = 171;                  % number of firms that draw
ebar    = 1.6379;               % productivity of type 2
su      = 0.1346;               % standard deviation of Home/Foreign productivity gap
tau     = 1.1263; 

om1     = [0.0778; 0.0150;  0.2808]; 
om2     = [0.1976; 0.2117;  0.1695]; 

pom     = [0.6312; 0.5759;  0.4427;... 
           0.6875; 0.5095;  0.5073;... 
           0.5697; 0.7446;  0.4649];
       
pom     = pom/sum(pom); 
      
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