function y = objective(x)

global gamma theta tau ebar N S

% Assigned Parameters
warning('off', 'MATLAB:nearlySingularMatrix')

gamma   = 10.4971;                    
theta   = 1.2351;
S       = 1000;                 % number of industries: must be even
L       = 1;                    % inelastic labor supply in Home
Ls      = 1;                    % inelastic labor supply in Foreign

lambda  = 0.0;                  % prob. all large guys collude
pom     = 1;                    % ok for now


% Calibrated

N       = x(1);                  % number of firms that draw
om1     = x(2); 
om2     = x(3); 
ebar    = x(4);                  % productivity of type 2
su      = x(5);                 % standard deviation of Home/Foreign productivity gap
tau     = x(6); 


     
comp_armington = 1;             % compute armington elasticity

[uu, wu] = qnwunif(7, 0, 1);
uu       = exp(norminv(uu, 0, su)); 

equilibrium
model_moments


loss(25)=100*loss(25);
loss=[loss;100*(-Arm-4)];
y = nanmean(abs(100*loss));

return;