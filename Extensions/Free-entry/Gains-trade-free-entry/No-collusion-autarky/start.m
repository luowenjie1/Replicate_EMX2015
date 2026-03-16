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

% Calibrated

N       = 168.0000;             % number of firms that draw
ebar    = 1.91481092395022;     % productivity of type 2

om1     = [0.0474556404890746; 0.0406677916576208; 0.293245448311502]; 
om2     = [0.180148410783646;  0.441073354951274;  0.421128903209272]; 

pom     = [0.586691719703357;  0.471324998862172;  0.621183930685627; 
           0.65345673615822;   0.869377396375917;  0.121756520349675; 
           0.40678208443039;   0.255111813402371;  0.259406018561665];
       
pom     = pom/sum(pom); 

%[uu, wu] = qnwunif(7, 0, 1);
%uu       = exp(norminv(uu, 0, su)); 

uu = 1; 
wu = 1; 

N = 187; 
%N = 168; 

equilibrium

