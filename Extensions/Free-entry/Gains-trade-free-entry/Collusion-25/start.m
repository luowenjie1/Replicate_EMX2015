clear all;
clc;

global gamma theta tau ebar N S

% Assigned Parameters

gamma   = 10.4971;                    
theta   = 1.2351;

S       = 1000;                 % number of industries: must be even
L       = 1;                    % inelastic labor supply in Home
Ls      = 1;                    % inelastic labor supply in Foreign

lambda  = 0.25;                 % prob. all large guys collude

% Calibrated

N       = 160;                  % number of firms that draw
ebar    = 1.7176;               % productivity of type 2
su      = 0.1737;               % standard deviation of Home/Foreign productivity gap
tau     = 1.1356; 

om1     = [0.1583; 0.0116; 0.3131]; 
om2     = [0.1535; 0.0587; 0.2596]; 

pom     = [0.4681; 0.3000; 0.4239; ...
           0.5701; 0.6364; 0.6787; ...
           0.4236; 0.7156; 0.6289];

pom     = pom/sum(pom); 

[uu, wu] = qnwunif(7, 0, 1);
uu       = exp(norminv(uu, 0, su)); 

equilibrium
