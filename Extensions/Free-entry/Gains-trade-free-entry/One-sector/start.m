clear all;
clc;

global gamma theta tau ebar N S

% Assigned Parameters

gamma   = 10.4971;                    
theta   = 1.2351;

S       = 1000;                 % number of industries: must be even
L       = 1;                    % inelastic labor supply in Home
Ls      = 1;                    % inelastic labor supply in Foreign

lambda  = 0;                    % prob. all large guys collude

% Calibrated

N       = 176;                  % number of firms that draw
om1     = 0.0972513900; 
om2     = 0.1734328984; 
ebar    = 1.838675631;          % productivity of type 2
su      = 0.1469919780;         % standard deviation of Home/Foreign productivity gap
tau     = 1.1339970969; 
             
pom     = 1; 

[uu, wu] = qnwunif(7, 0, 1);
uu       = exp(norminv(uu, 0, su)); 


equilibrium
