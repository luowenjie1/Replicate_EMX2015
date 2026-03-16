rng(0);   % reset random number generator

% Idiosyncratic firm productivity

[ux,uy] = gumbel_copula(N*S,0); %% copula function, independent since kendall's tau=0 

a  = (1-ux).^(-1/xi_x);
as = (1-uy).^(-1/xi_x);

clear ux uy; 

nn = geornd(zeta,S,1); %% sample from geometric distribution with parameter zeta


phi = zeros(N, S); 

for i = 1:S
    
    temp = zeros(N,1); 
    temp(1:nn(i)) = 1; 
    temp = temp(randperm(N)); 
    phi(:,i) = temp; 
    
end

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
as = as.*phi + eps.*(1-phi);


a  = sort(a, 1,'descend');
as = sort(as,1,'descend');

% MAKE COUNTRIES PERFECTLY SYMMETRIC: Drop 1/2 observations, switch around

a  = [a, as];
as = [a(:, S+1:end), a(:,1:S)];


S  = S*2;

a  = sort(a, 1, 'descend');
as = sort(as, 1, 'descend');


% initial guess: economy with p = gamma/(gamma-1)*W/a

tariff=[ones(N,S);1-repmat(betarnd(tariff_a,tariff_b,1,S),N,1)]; %%(1-tariff)

aa = [a;as/(1+tau)]; 

pp  = gamma/(gamma-1)./(aa.*tariff);


Pj  = (mean(2*pp.^(1-gamma))).^(1/(1-gamma));
P   = (mean(Pj.^(1-theta))).^(1/(1-theta));

phi = ones(2*N, S);

A  = mean(mean(2*phi.*aa.^(-1).*(pp./repmat(Pj,2*N,1)).^(-gamma).*repmat(Pj/P, 2*N,1).^(-theta))).^(-1);
     
Y     = A*L;
Yold  = Y;

% iterate until convergence on firm decision rules and exit decisions

omegaold  =  1/N*(pp./repmat(Pj,2*N,1)).^(1-gamma); 
yyold      = (pp./repmat(Pj,2*N,1)).^(-gamma).*(repmat(Pj,2*N,1)./P).^(-theta).*Y;


for it = 1 : 200
   
Y      = Yold;  
omega  = omegaold;
yy     = yyold;

ee = (1/gamma*(1-omega) + 1/theta*omega).^(-1); 
pp  = ee./(ee-1)./(aa.*tariff);

profit = (tariff.*pp - 1./aa).*yy - [FD*ones(N,S); FX*ones(N,S)];

if it < 50    
  
phi  = profit >= 0;   % don't iterate on phi further, discreteness gives jumps, no convergence

end

Pj              = mean(2*phi.*pp.^(1-gamma)).^(1/(1-gamma));
Pj(isinf(Pj))   = 10^16;

P               = mean(Pj.^(1-theta)).^(1/(1-theta));

omega = 1/N*bsxfun(@times, pp.^(1-gamma).*phi, Pj.^(gamma-1));

yy = bsxfun(@times, pp.^(-gamma), Pj.^(gamma-theta).*P.^theta.*Y);

omega = real(omega);
yy    = real(yy); 

A  = mean(mean(2*phi.*aa.^(-1).*yy/Y)).^(-1);   

fcost = phi.*[FD*ones(N,S); FX*ones(N,S)];
Y     = A*(L-2*mean(fcost(:)));

error = max([norm(omega - omegaold); norm(Y-Yold)]);


fprintf('%4i %6.2e \n',[it, error]);


if error < 1e-7, break, end
    
Yold         = 1/2*Y      +  1/2*Yold;
omegaold     = 1/10*omega +  9/10*omegaold;
yyold        = 1/10*yy    +  9/10*yyold;

end

yy  = yy.*(phi>0);
omega = omega.*(phi>0);

ll = yy./aa;

if norm(omega-omegaold)>1e-3
    disp('Warning: omega has not converged')
end

mu  = ee./(ee-1)./tariff;  

% Naive Trade Elasticity (fix measures of producers)

omegaH = omega(1:N,:); 
omegaF = omega(N+1:2*N,:); 
pH     = pp(1:N,:); 
pF     = pp(N+1:2*N,:);    % this is (1+tau)*pF
yH     = yy(1:N,:); 
yF     = yy(N+1:2*N,:); 
phiH   = phi(1:N,:); 
phiF   = phi(N+1:2*N,:);
muH    = mu(1:N,:); 
muF    = mu(N+1:2*N,:); 
lH     = ll(1:N,:); 
lF     = ll(N+1:2*N,:);


% domestic and foreign markups
mudom = (mean(mean(muH.^(-1).*pH.*yH))/mean(mean(pH.*yH)))^(-1);
mufor = (mean(mean((1+tau).*muF.^(-1).*pF.*yF))/mean(mean((1+tau)*pF.*yF)))^(-1);
domshare = mean(mean(pH.*yH))/(P*Y);

if domshare>0.999,
    
    mucheck = mudom; %% autarky, mufor not defined
    
else

    mucheck = (mudom^(-1)*domshare + mufor^(-1)*(1-domshare))^(-1); %% should be same as aggregate markup

end
    
muagg   = P*A;


lambdaj = sum(omegaH);
sj      = (Pj./P).^(1-theta);
lambda  = mean(lambdaj.*sj);
aweight = mean(sj.*lambdaj/lambda.*(1-lambdaj)/(1-lambda));
naiveArm     = gamma*aweight + theta*(1-aweight);

if abs(lambda-1) < 1e-3   % don't divide by 0
  
    aweight = 1;
    naiveArm = gamma;
    
end

relmarkups  = mean(mean(phiF.*(pF./repmat(Pj,N,1)).^(1-gamma).*(repmat(Pj,N,1)/P).^(1-theta)))/...
              mean(mean(phiH.*(pH./repmat(Pj,N,1)).^(1-gamma).*(repmat(Pj,N,1)/P).^(1-theta)));

impshare    = sum(omegaF)';
import_pen  = impshare; 

% use the fact that everything is symmetric (a(:,j) = as(:,S/2+j), as(:,j)
% = a(:,S/2+j), to back out "export share" of industry j from 
% = import share of industry S/2 + j

pHs     = [pF(:,S/2+1:S), pF(:,1:S/2)];
yHs     = [yF(:,S/2+1:S), yF(:,1:S/2)];
phiHs   = [phiF(:,S/2+1:S), phiF(:,1:S/2)];



omegaHs = [omegaF(:,S/2+1:S), omegaF(:,1:S/2)];

expshare   = sum(omegaHs)';

ttrade     = import_pen + expshare;
intraindex  = mean(sj(ttrade>0)'.*abs(impshare(ttrade>0)-expshare(ttrade>0))./(impshare(ttrade>0) + expshare(ttrade>0)))./mean(sj(ttrade>0)');
     
impshare   = sum(pF(:).*yF(:))/sum(pH(:).*yH(:) + pF(:).*yF(:));

agg_fexporters = sum(omegaH(:)>0 & omegaF(:)>0)/sum(omegaH(:)>0);

% Correct Trade Elasticity (trade costs affect markups and measures of
% producers)

if comp_trade_elasticity
    
parameters   = [gamma, theta, N, S, FD, FX ];

step   = 1.005;  
taunew = (1+tau)*step - 1;
aanew  = [a; as/(1+taunew)]; 

sigma = trade_elasticity(aanew, Y, P, omega, yy, parameters , step, impshare , tariff);

%Arm = log(impsharenew./(1-impsharenew)/(impshare/(1-impshare)))/log(step);

else
    
sigma = 0;
    
end

% just use the same formulas as before but with
% different prices

ppe             = 1./(aa.*tariff);  % we now assume all firms charge a constant markup, here normalized to 1 since mean irrelevant
Pje             = mean(2*phi.*ppe.^(1-gamma)).^(1/(1-gamma));
Pje(isinf(Pje)) = 10^16;
Pe              = mean(Pje.^(1-theta)).^(1/(1-theta));

yye   = bsxfun(@times, ppe.^(-gamma), Pje.^(gamma-theta).*Pe.^theta); 
Aeff  = mean(mean(2*phi.*aa.^(-1).*yye)).^(-1);                         


sshare = (Pj./P).^(1-theta); sshare = sshare';

Lj   = sum(lH+lF); 
PjYj = sum(pH.*yH + pF.*yF); 
muj  = PjYj'./Lj';

fprintf('\n');
display('Equilibrium Implications');
fprintf('\n');
fprintf('Y                 = %7.3f \n', Y);  
fprintf('A                 = %7.3f \n', A);
fprintf('A loss, *100      = %7.3f \n', log(Aeff/A)*100);  
fprintf('mean markup       = %7.3f \n', nanmean(mu(phi>0).^(-1)).^(-1));
fprintf('mean sect. markup = %7.3f \n', nanmean(muj.^(-1)).^(-1));
fprintf('aggregate markup  = %7.3f \n', P*A);
fprintf('import share      = %7.3f \n', impshare);
fprintf('trade elasticity  = %7.3f \n', sigma);
fprintf('\n');
fprintf('domestic markup   = %7.3f \n', mudom);
fprintf('import   markup   = %7.3f \n', mufor);
fprintf('\n');



clear a as
S = S/2;




