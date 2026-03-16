rng(0);              % reset random number generator

% Idiosyncratic firm productivity

[ux,uy] = gumbel_copula(N*S,0); %% copula function, independent since Kendall tau=0

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

a  = sort(a,  1, 'descend');
as = sort(as, 1, 'descend');

% initial guess: economy with p = gamma/(gamma-1)*W/a

Ws  = 1;     % guess wage in foreign

aa  = [a;            1/(1+tau)*as*Abars/Ws]; 
aas = [as*Abars/Ws;    1/(1+tau)*a       ];

mu  = gamma/(gamma-1); 
mus = gamma/(gamma-1); 

Pj  = (1/N*sum(mu.^(1-gamma).*aa.^(gamma-1))).^(1/(1-gamma));
Pjs = (1/N*sum(mus.^(1-gamma).*aas.^(gamma-1))).^(1/(1-gamma));
P   = (mean(Pj.^(1-theta))).^(1/(1-theta));
Ps  = (mean(Pjs.^(1-theta))).^(1/(1-theta));

eta = 1; % from balanced trade equation

phi  = ones(2*N, S);
phis = ones(2*N, S); 

A    = mean(1/N*sum(phi(1:N,:).*mu.^(-gamma).*aa(1:N,:).^(gamma-1)).*Pj.^(gamma-theta)*P^theta + ...
          eta/N*sum(phis(N+1:2*N,:).*mus.^(-gamma).*aas(N+1:2*N,:).^(gamma-1)).*Pjs.^(gamma-theta)*P^theta).^(-1);

As   = (1/Ws*mean(1/N*sum(phis(1:N,:).*mus.^(-gamma).*aas(1:N,:).^(gamma-1)).*Pjs.^(gamma-theta)*Ps^theta + ...
          eta^(-1)/N*sum(phi(N+1:2*N,:).*mu.^(-gamma).*aa(N+1:2*N,:).^(gamma-1)).*Pj.^(gamma-theta)*Ps^theta)).^(-1);
             
Y     = A*L;
Ys    = As*Ls;

Yold  = Y;
Ysold = Ys; 
Wsold = Ws; 

% iterate until convergence on firm decision rules and exit decisions

omegaold  =  1/N*phi.*(mu./aa./repmat(Pj,2*N,1)).^(1-gamma); 
omegasold =  1/N*phis.*(mus./aas./repmat(Pjs,2*N,1)).^(1-gamma); 

yyold      = phi.*(mu./aa).^(-gamma).*repmat(Pj,2*N,1).^(gamma-theta).*P.^(theta).*Y; % output
yysold     = phis.*(mus./aas).^(-gamma).*repmat(Pjs,2*N,1).^(gamma-theta).*Ps.^(theta).*Ys; % output

for it = 1 : 500
   
Y      = Yold;  
Ys     = Ysold; 
Ws     = Wsold; 

omega  = omegaold;
omegas = omegasold; 

yy     = yyold;
yys    = yysold; 

aa  = [a;           1/(1+tau)*as*Abars/Ws]; 
aas = [as*Abars/Ws;    1/(1+tau)*a       ];

ee  = (1/gamma*(1-omega) + 1/theta*omega).^(-1); 
ees = (1/gamma*(1-omegas) + 1/theta*omegas).^(-1); 

mu  = ee./(ee-1);
mus = ees./(ees-1);

pp  = mu./aa;
pps = mus./aas;

profit  = (pp - 1./aa).*yy - [FD*ones(N,S); FX*ones(N,S)];
profits = (pps - 1./aas).*yys - [FDs*ones(N,S); FX*ones(N,S)];

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

eta    = mean(1/N*sum(phi(N+1:2*N,:).*mu(N+1:2*N,:).^(1-gamma).*aa(N+1:2*N,:).^(gamma-1)).*Pj.^(gamma-theta))/...
         mean(1/N*sum(phis(N+1:2*N,:).*mus(N+1:2*N,:).^(1-gamma).*aas(N+1:2*N,:).^(gamma-1)).*Pjs.^(gamma-theta));
     
A      = mean(1/N*sum(phi(1:N,:).*mu(1:N,:).^(-gamma).*aa(1:N,:).^(gamma-1)).*Pj.^(gamma-theta)*P^theta + ...
         eta/N*sum(phis(N+1:2*N,:).*mus(N+1:2*N,:).^(-gamma).*aas(N+1:2*N,:).^(gamma-1)).*Pjs.^(gamma-theta)*P^theta).^(-1);

As     = Ws*mean(1/N*sum(phis(1:N,:).*mus(1:N,:).^(-gamma).*aas(1:N,:).^(gamma-1)).*Pjs.^(gamma-theta)*Ps^theta + ...
          eta^(-1)/N*sum(phi(N+1:2*N,:).*mu(N+1:2*N,:).^(-gamma).*aa(N+1:2*N,:).^(gamma-1)).*Pj.^(gamma-theta)*Ps^theta).^(-1);
  
Ws     = Abars*(mean(1/N*sum(phis(N+1:2*N,:).*mus(N+1:2*N,:).^(1-gamma).*a.^(gamma-1)).*Pjs.^(gamma-theta)*Ps^theta*Ys)/...
           mean(1/N*sum(phi(N+1:2*N,:).*mu(N+1:2*N,:).^(1-gamma).*as.^(gamma-1)).*Pj.^(gamma-theta)*P^theta*Y) ).^(1/(1-gamma));       

fcost  = phi.*[FD*ones(N,S); FX*ones(N,S)];
Y      = A*(L-2*mean(fcost(:)));

fcosts = phis.*[FDs*ones(N,S); FX*ones(N,S)];
Ys     = As*(Ls-2*mean(fcosts(:)));

error  = max([norm(omega - omegaold); norm(omegas - omegasold); norm([Y-Yold; Ys - Ysold; Ws - Wsold])]);

fprintf('%4i %6.2e \n', [it, error]);

if error < 1e-7, break, end
    
Yold         = 1/2*Y    +  1/2*Yold;
Ysold        = 1/2*Ys   +  1/2*Ysold; 
Wsold        = 1/2*Ws   +  1/2*Wsold; 

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


omegaH  = omega(1:N,:); 
omegaF  = omega(N+1:2*N,:); 
omegaHs = omegas(N+1:2*N,:); 
omegaFs = omegas(1:N,:);

pH      = pp(1:N,:); 
pF      = pp(N+1:2*N,:);    % this is (1+tau)*pF
pHs     = pps(N+1:2*N,:);
pFs     = pps(1:N,:);

yH      = yy(1:N,:); 
yF      = yy(N+1:2*N,:); 
yHs     = yys(N+1:2*N,:);
yFs     = yys(1:N,:);

phiH    = phi(1:N,:); 
phiF    = phi(N+1:2*N,:);
phiHs   = phis(N+1:2*N,:);
phiFs   = phis(1:N,:); 

muH     = mu(1:N,:); 
muF     = mu(N+1:2*N,:); 
muHs    = mus(N+1:2*N,:);
muFs    = mus(1:N,:);

lH      = yy(1:N,:)./aa(1:N,:); 
lF      = yy(N+1:2*N,:)./aa(N+1:2*N,:)/Ws;   % notice it includes (1+tau) since in a, but must take out Ws
lHs     = yys(N+1:2*N,:)./aas(N+1:2*N,:);
lFs     = yys(1:N,:)./aas(1:N,:)/Ws;

impshare   = sum(pF(:).*yF(:))/sum(pH(:).*yH(:) + pF(:).*yF(:));
impshares  = sum(pHs(:).*yHs(:))/sum(pFs(:).*yFs(:) + pHs(:).*yHs(:));
impsharea  = sum(pF(:).*yF(:) + pHs(:).*yHs(:))/sum(pH(:).*yH(:) + pF(:).*yF(:) + pFs(:).*yFs(:) + pHs(:).*yHs(:) );

profitH = profit(1:N,:).*phiH; 
profitF = profit(N+1:2*N,:).*phiF;  
profitFs = profits(1:N,:).*phiFs;
profitHs = profits(N+1:2*N,:).*phiHs; 


sj  = (Pj./P).^(1-theta);
sjs = (Pjs./Ps).^(1-theta);

lambdaj  = sum(omegaH);
lambdajs = sum(omegaFs); 

lambda   = mean(lambdaj.*sj);
lambdas  = mean(lambdajs.*sjs);

aweight  = mean(sj.*lambdaj/lambda.*(1-lambdaj)/(1-lambda));
aweights = mean(sjs.*lambdajs/lambdas.*(1-lambdajs)/(1-lambdas));


if comp_trade_elasticity

parameters   = [gamma, theta, N, S, FD, FX];
    
step   = 1.005; 
taunew = (1+tau)*step - 1;

aanew  = [a;              1/(1+taunew)*as*Abars/Ws ]; 
aasnew = [as*Abars/Ws;    1/(1+taunew)*a           ]; 

[sigma, sigmas] = trade_elasticity(aanew, aasnew, Y, Ys, P, Ps, omega, omegas, yy, yys, parameters, step, impshare, impshares); 

%Arm  = log(impsharenew./(1-impsharenew)/(impshare/(1-impshare)))/log(step);
%Arms = log(impsharesnew./(1-impsharesnew)/(impshares/(1-impshares)))/log(step);

else
    
sigma = 0;
sigmas = 0;
    
end

Aje = mean(2*phi.*aa.^(gamma-1)).^(1/(gamma-1)); 
Ae  = mean(Aje.^(theta-1)).^(1/(theta-1)); 

Ajse = mean(2*phis.*aas.^(gamma-1)).^(1/(gamma-1)); 
Ase  = Ws*mean(Ajse.^(theta-1)).^(1/(theta-1)); 



      

Lj   = sum(lH+lF); 
PjYj = sum(pH.*yH + pF.*yF); 
muj  = PjYj'./Lj';

Ljs   = sum(lHs+lFs); 
PjYjs = sum(pHs.*yHs + pFs.*yFs); 
mujs  = PjYjs'./(Ws*Ljs)';
      
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

clear a as
S=S/2;

agg_fexporters  = sum(omegaH(:)>0 & omegaF(:)>0)/sum(omegaH(:)>0); 
agg_fexporterss = sum(omegaHs(:)>0 & omegaFs(:)>0)/sum(omegaFs(:)>0); 