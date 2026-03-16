clear all;
clc;

global gamma theta ebar N S

% Assigned Parameters

gamma   = 10.4971;                    
theta   = 1.2351;

S       = 1000;                 % number of industries: must be even
L       = 1;                    % inelastic labor supply in Home
Ls      = 1;                    % inelastic labor supply in Foreign

lambda  = 0;                    % prob. all large guys collude

N       = 176;                  % number of firms that draw
om1     = 0.0972513900; 
om2     = 0.1734328984; 
ebar    = 1.838675631;          % productivity of type 2
su      = 0.1469919780;         % standard deviation of Home/Foreign productivity gap
tau     = 1.1339970969; 
             
pom     = 1; 

uu = 1; 
wu = 1; 

N = 191; 

equilibrium
