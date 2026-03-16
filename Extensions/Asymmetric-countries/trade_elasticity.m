function [sigma, sigmas] = trade_elasticity(aa, aas, Y, Ys, P, Ps,  omega, omegas, yy, yys,  parameters , step, impshare, impshares)

gamma = parameters(1);
theta = parameters(2);
N     = parameters(3);
S     = parameters(4);
FD    = parameters(5);
FX    = parameters(6);
 
omegaold  = omega;
omegasold = omegas; 

yyold    = yy; 
yysold   = yys; 

for it = 1 : 200
   
omega  = omegaold;
omegas = omegasold; 

yy     = yyold;
yys    = yysold; 

ee  = (1/gamma*(1-omega) + 1/theta*omega).^(-1); 
ees = (1/gamma*(1-omegas) + 1/theta*omegas).^(-1); 

mu  = ee./(ee-1);
mus = ees./(ees-1);

pp  = mu./aa;
pps = mus./aas;

profit  = (pp - 1./aa).*yy - [FD*ones(N,S); FX*ones(N,S)];
profits = (pps - 1./aas).*yys - [FD*ones(N,S); FX*ones(N,S)];

if it < 50    
  
phi  = profit >= 0;   % don't iterate on phi further, discreteness gives jumps, no convergence
phis = profits >= 0; 

end

Pj  = (1/N*sum(phi.*mu.^(1-gamma).*aa.^(gamma-1))).^(1/(1-gamma));
Pj(isinf(Pj))   = 10^16;

Pjs = (1/N*sum(phis.*mus.^(1-gamma).*aas.^(gamma-1))).^(1/(1-gamma));
Pjs(isinf(Pjs))   = 10^16;

omega  = 1/N*bsxfun(@times, pp.^(1-gamma).*phi, Pj.^(gamma-1));
omegas = 1/N*bsxfun(@times, pps.^(1-gamma).*phis, Pjs.^(gamma-1));

yy     = bsxfun(@times, pp.^(-gamma),  Pj.^(gamma-theta).*P.^theta.*Y);
yys    = bsxfun(@times, pps.^(-gamma), Pjs.^(gamma-theta).*Ps.^theta.*Ys);

omega  = real(omega);
omegas = real(omegas);

yy     = real(yy); 
yys    = real(yys); 

error  = max([norm(omega - omegaold); norm(omegas-omegasold)]);

if error < 1e-7, break, end
    
omegaold   = 1/10*omega  +  9/10*omegaold;
omegasold  = 1/10*omegas +  9/10*omegasold; 

yyold      = 1/10*yy    +  9/10*yyold;
yysold     = 1/10*yys   +  9/10*yysold; 

end

yy    = yy.*(phi>0);
yys   = yys.*(phis>0); 

pH      = pp(1:N,:); 
pF      = pp(N+1:2*N,:);    % this is (1+tau)*pF
pHs     = pps(N+1:2*N,:);
pFs     = pps(1:N,:);

yH      = yy(1:N,:); 
yF      = yy(N+1:2*N,:); 
yHs     = yys(N+1:2*N,:);
yFs     = yys(1:N,:);

impsharenew   = sum(pF(:).*yF(:))/sum(pH(:).*yH(:) + pF(:).*yF(:));
impsharesnew  = sum(pHs(:).*yHs(:))/sum(pFs(:).*yFs(:) + pHs(:).*yHs(:));


sigma  = -log(impsharenew./(1-impsharenew)/(impshare/(1-impshare)))/log(step);
sigmas = -log(impsharesnew./(1-impsharesnew)/(impshares/(1-impshares)))/log(step);