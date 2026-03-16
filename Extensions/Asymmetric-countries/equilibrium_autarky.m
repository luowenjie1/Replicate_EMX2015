rng(0);              % reset random number generator

% Idiosyncratic firm productivity

%seed = rng

[ux,uy] = gumbel_copula(N*S,0); %% copula function, taurhoi=0

a  = (1-ux).^(-1/xi_x);
as = (1-uy).^(-1/xi_x);

clear ux uy; 


nn = geornd(zeta,S,1); %% sample from geometric distribution with parameter zeta

phi  = zeros(N, S); 
phis = zeros(N,S);

for i = 1:S
    
    temp = zeros(N,1); 
    temp(1:nn(i)) = 1; 
    temp = temp(randperm(N)); 
    phi(:,i) = temp; 
    
end

phis = phi;

a   = reshape(a, N, S);
as  = reshape(as, N, S); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sector productivity 

[ux,uy] = gumbel_copula(S,rho); %% copula function

z  = (1-ux).^(-1/xi_z);
zs = (1-uy).^(-1/xi_z);

clear ux uy; 

a  = bsxfun(@times, a, z'); 
as = bsxfun(@times, as, zs'); 

a  = a.*phi  + eps.*(1-phi);
as = as.*phis + eps.*(1-phis);


a  = sort(a, 1,'descend');
as = sort(as,1,'descend');

% MAKE COUNTRIES PERFECTLY SYMMETRIC: Drop 1/2 observations, switch around

a  = [a, as];
as = [a(:, S+1:end), a(:,1:S)];

S  = S*2;

a  = sort(a, 1, 'descend');
as = sort(as, 1, 'descend');

% initial guess: economy with p = gamma/(gamma-1)*W/a

Ws  = 1;     % guess wage in foreign (it is no longer pinned down, since autarky but check)

aa  = a; 
aas = as*Abars/Ws; 

mu  = gamma/(gamma-1); 
mus = gamma/(gamma-1); 

Pj  = (1/N*sum(mu.^(1-gamma).*aa.^(gamma-1))).^(1/(1-gamma));
Pjs = (1/N*sum(mus.^(1-gamma).*aas.^(gamma-1))).^(1/(1-gamma));
P   = (mean(Pj.^(1-theta))).^(1/(1-theta));
Ps  = (mean(Pjs.^(1-theta))).^(1/(1-theta));

phi  = ones(N, S);
phis = ones(N, S); 

A    = mean(1/N*sum(phi(1:N,:).*mu.^(-gamma).*aa(1:N,:).^(gamma-1)).*Pj.^(gamma-theta)*P^theta).^(-1);

As   = (1/Ws*mean(1/N*sum(phis(1:N,:).*mus.^(-gamma).*aas(1:N,:).^(gamma-1)).*Pjs.^(gamma-theta)*Ps^theta)).^(-1);
             
Y     = A*L;
Ys    = As*Ls;

Yold  = Y;
Ysold = Ys; 
Wsold = Ws; 

% iterate until convergence on firm decision rules and exit decisions

omegaold  =  1/N*phi.*(mu./aa./repmat(Pj,N,1)).^(1-gamma); 
omegasold =  1/N*phis.*(mus./aas./repmat(Pjs,N,1)).^(1-gamma); 

yyold      = phi.*(mu./aa).^(-gamma).*repmat(Pj,N,1).^(gamma-theta).*P.^(theta).*Y; % output
yysold     = phis.*(mus./aas).^(-gamma).*repmat(Pjs,N,1).^(gamma-theta).*Ps.^(theta).*Ys; % output


for it = 1 : 200
   
Y      = Yold;  
Ys     = Ysold; 
Ws     = Wsold; 

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

profit  = (pp - 1./aa).*yy - FD*ones(N,S);
profits = (pps - 1./aas).*yys - FD*ones(N,S);

if it < 50    
  
phi  = profit >= 0;   % don't iterate on phi further, discreteness gives jumps, no convergence
phis = profits >= 0; 

end

Pj     = (1/N*sum(phi.*mu.^(1-gamma).*aa.^(gamma-1))).^(1/(1-gamma));
Pj(isinf(Pj))   = 10^16;

Pjs    = (1/N*sum(phis.*mus.^(1-gamma).*aas.^(gamma-1))).^(1/(1-gamma));
Pjs(isinf(Pjs))   = 10^16;

P      = (mean(Pj.^(1-theta))).^(1/(1-theta));
Ps     = (mean(Pjs.^(1-theta))).^(1/(1-theta));

omega  = 1/N*bsxfun(@times, pp.^(1-gamma).*phi, Pj.^(gamma-1));
omegas = 1/N*bsxfun(@times, pps.^(1-gamma).*phis, Pjs.^(gamma-1));

yy     = bsxfun(@times, pp.^(-gamma),  Pj.^(gamma-theta).*P.^theta.*Y);
yys    = bsxfun(@times, pps.^(-gamma), Pjs.^(gamma-theta).*Ps.^theta.*Ys);

omega  = real(omega);
omegas = real(omegas);

yy     = real(yy); 
yys    = real(yys); 

A      = mean(1/N*sum(phi(1:N,:).*mu(1:N,:).^(-gamma).*aa(1:N,:).^(gamma-1)).*Pj.^(gamma-theta)*P^theta).^(-1);

As     = Ws*mean(1/N*sum(phis(1:N,:).*mus(1:N,:).^(-gamma).*aas(1:N,:).^(gamma-1)).*Pjs.^(gamma-theta)*Ps^theta).^(-1);
  
fcost  = phi.*FD.*ones(N,S);
Y      = A*(L-2*mean(fcost(:)));

fcosts = phis.*FD.*ones(N,S);
Ys     = As*(Ls-2*mean(fcosts(:)));

error  = max([norm(omega - omegaold); norm(omegas - omegasold); norm([Y-Yold; Ys - Ysold])]);


fprintf('%4i %6.2e \n', [it, error]);


if error < 1e-7, break, end
    
Yold         = 1/2*Y    +  1/2*Yold;
Ysold        = 1/2*Ys   +  1/2*Ysold; 

omegaold     = 1/10*omega  +  9/10*omegaold;
omegasold    = 1/10*omegas +  9/10*omegasold; 

yyold        = 1/10*yy    +  9/10*yyold;
yysold       = 1/10*yys   +  9/10*yysold; 

end

yy     = yy.*(phi>0);
omega  = omega.*(phi>0);

yys    = yys.*(phis>0); 
omegas = omegas.*(phis>0);  


if norm(omega-omegaold) > 1e-3 || norm(omegas-omegasold) > 1e-3
    disp('Warning: omega has not converged')
end



% losses from misallocation

Pje               = (1/N*sum(phi.*aa.^(gamma-1))).^(1/(1-gamma));
Pje(isinf(Pje))   = 10^16;

Pjse              = (1/N*sum(phis.*aas.^(gamma-1))).^(1/(1-gamma));
Pjse(isinf(Pjse)) = 10^16;

Pe   = (mean(Pje.^(1-theta))).^(1/(1-theta));
Pse  = (mean(Pjse.^(1-theta))).^(1/(1-theta));


Ae     =  mean(1/N*sum(phi(1:N,:).*aa(1:N,:).^(gamma-1)).*Pje.^(gamma-theta)*Pe^theta).^(-1);

Ase    =  Ws*mean(1/N*sum(phis(1:N,:).*aas(1:N,:).^(gamma-1)).*Pjse.^(gamma-theta)*Pse^theta).^(-1);
 
fprintf('\n');
display('Equilibrium Implications');
fprintf('\n');
fprintf('Y                  = %7.3f \n', Y);  
fprintf('Ys                 = %7.3f \n', Ys);  
fprintf('\n');
fprintf('A                  = %7.3f \n', A);
fprintf('As                 = %7.3f \n', As);
fprintf('\n');
fprintf('A loss, *100       = %7.3f \n', log(Ae/A)*100);  
fprintf('As loss, *100      = %7.3f \n', log(Ase/As)*100);  
fprintf('\n');

fprintf('import share H     = %7.3f \n', 0);
fprintf('import share F     = %7.3f \n', 0);
fprintf('\n');

clear a as 
S = S/2;