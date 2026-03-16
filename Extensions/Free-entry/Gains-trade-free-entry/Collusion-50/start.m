clear all;
clc;

global gamma theta tau ebar N S

% Assigned Parameters

gamma   = 10.4971;                    
theta   = 1.2351;

S       = 1000;                 % number of industries: must be even
L       = 1;                    % inelastic labor supply in Home
Ls      = 1;                    % inelastic labor supply in Foreign

lambda  = 0.50;                 % prob. all large guys collude

% Calibrated

N       = 110;                  % number of firms that draw
ebar    = 1.4866;               % productivity of type 2
su      = 0.1955;               % standard deviation of Home/Foreign productivity gap
tau     = 1.1331; 

om1     = [0.1216; 0.0283; 0.3388]; 
om2     = [0.3271; 0.1262; 0.2030]; 

pom     = [0.4775; 0.8366; 0.1584;... 
           0.2276; 0.2994; 0.7672;... 
           0.9008; 0.7324; 0.0385];
       
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