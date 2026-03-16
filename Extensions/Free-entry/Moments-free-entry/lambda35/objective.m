function y = objective(x)

global gamma theta tau ebar N S

% Assigned Parameters
warning('off', 'MATLAB:nearlySingularMatrix')

gamma   = 10.4971;                    
theta   = 1.2351;

S       = 1000;                 % number of industries: must be even
L       = 1;                    % inelastic labor supply in Home
Ls      = 1;                    % inelastic labor supply in Foreign

lambda  = 0.35;                  % prob. all large guys collude
pom     = 1;                    % ok for now


% Calibrated

N       = x(1);                  % number of firms that draw
ebar    = x(2);                  % productivity of type 2
su      = x(3);                 % standard deviation of Home/Foreign productivity gap
tau     = x(4); 

om1     = [x(5); x(6); x(7)]; 
om2     = [x(8); x(9); x(10)]; 

pom     = [x(11);x(12);x(13);
           x(14);x(15);x(16);
           x(17);x(18);x(19)]; 

save It_x x


pom=pom./sum(pom);

comp_armington = 1;             % compute armington elasticity

[uu, wu] = qnwunif(7, 0, 1);
uu       = exp(norminv(uu, 0, su)); 

equilibrium
display(x)
model_moments


loss(25)=100*loss(25);
loss=[loss;100*(-Arm-4)];
y = nanmean(abs(100*loss));
display(x)
display(y)
return;