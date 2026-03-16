%rng(0);   % reset random number generator
rand('state', 0) 
% For each industry type (Home/Foreign productivity), draw n1(i), n2(i),
% n1*(i), n2*(i)

if 0

xi  = [1; 0];  
xis = [1; 0]; 
    
typei = gridmake(om1, om2, uu, xi, xis); 
probi = gridmake(pom, wu, [lambda; 1-lambda], [lambda; 1-lambda]);
probi = probi(:,1).*probi(:,2).*probi(:,3).*probi(:,4); 

else

    xi  = 0; 
    xis = 0; 
    
typei = gridmake(om1, om2, uu, xi, xis); 
probi = gridmake(pom, wu, 1, 1);
probi = probi(:,1).*probi(:,2).*probi(:,3).*probi(:,4); 

end    

K = size(typei,1); 

n1all  = zeros(S, K); 
n2all  = n1all;

for  itt = 1 : K
    
n            = binornd(N, typei(itt,1), S, 1) + 1;  % there is at least one firm in an industry
n2all(:,itt) = binornd(n, typei(itt,2));
n1all(:,itt) = n - n2all(:,itt); 

end

clear n

mu1Hall = zeros(S, K);
mu2Hall = zeros(S, K); 
mu1Fall = zeros(S, K);
mu2Fall = zeros(S, K); 
Piall   = zeros(S, K); 

fprintf('\n')
fprintf('Solve each industrys equilibrium \n')
fprintf('\n')

for itt = 1 : K
        
n1  = n1all(:,itt); 
n2  = n2all(:,itt); 
    
ui  = typei(itt,3); 
xi  = typei(itt,4);   % can only collude if n2>1
xis = typei(itt,5); 
    
n1s = [n1(S/2+1 : S); n1(1 : S/2)]; 
n2s = [n2(S/2+1 : S); n2(1 : S/2)];        

xi   = xi*(n2>1);       % collude dummy
xis  = xis*(n2s>1);     % 4 cases here, eventually list them all

% compute industry equilibrium for a given u(i), xi(i), xis(i)

mu1H = gamma/(gamma-1)*ones(S,1).*(n1>0) + (n1==0); 
mu2H = gamma/(gamma-1)*ones(S,1).*(n2>0) + (n2==0); 
mu1F = gamma/(gamma-1)*ones(S,1).*(n1s>0) + (n1s==0); 
mu2F = gamma/(gamma-1)*ones(S,1).*(n2s>0) + (n2s==0);

Pi   = (n1.*(mu1H/ui).^(1-gamma)   + n2.*(mu2H/ui/ebar).^(1-gamma) + ...
        n1s.*(tau*mu1F).^(1-gamma) + n2s.*(tau*mu2F/ebar)).^(1/(1-gamma));

xold = [mu1H, mu2H, mu1F, mu2F, Pi];     

err  = 1000; 

% A few function iterations

for it = 1: 25
   
mu1H = xold(:,1); 
mu2H = xold(:,2); 
mu1F = xold(:,3); 
mu2F = xold(:,4); 
Pi   = xold(:,5); 

    mu1H = (1 - 1/gamma + (mu1H./Pi).^(1-gamma).*(1/gamma-1/theta).*ui^(gamma-1)).^(-1).*(n1>0) + (n1==0);   
    mu2H = (max(1 - 1/gamma + (mu2H./Pi).^(1-gamma).*(1/gamma-1/theta).*(ui*ebar)^(gamma-1).*(n2.*xi + (1-xi)), 1-1/theta)).^(-1).*(n2>0) + (n2==0);
    mu1F = (1 - 1/gamma + (mu1F./Pi).^(1-gamma).*(1/gamma-1/theta).*(1/tau)^(gamma-1)).^(-1).*(n1s>0) + (n1s==0);   
    mu2F = (max(1 - 1/gamma + (mu2F./Pi).^(1-gamma).*(1/gamma-1/theta).*(ebar/tau)^(gamma-1).*(n2s.*xis + (1-xis)), 1-1/theta)).^(-1).*(n2s>0) + (n2s==0);
    Pi   = (n1.*(mu1H/ui).^(1-gamma) + n2.*(mu2H/ui/ebar).^(1-gamma) + ...
            n1s.*(tau*mu1F).^(1-gamma) + n2s.*(tau*mu2F/ebar).^(1-gamma)).^(1/(1-gamma));

xnew =  [mu1H, mu2H, mu1F, mu2F, Pi];     

%fprintf('%4i %6.2e \n',[it, norm(xnew - xold)]);   %report difference

err  = norm(xnew - xold);

xold = 1/20*xnew + 19/20*xold; 

end
   
x = [1./mu1H, 1./mu2H, 1./mu1F, 1./mu2F, Pi.^(gamma-1)];

% Newton

it = 0; 

while err > 1e-7 && it < 250
    
    it = it + 1; 
    
[y, yd] = solveequil(x, tau, ui, xi, xis, n1, n2, n1s, n2s, gamma, theta, ebar);

dx = - arrayinv(y, yd); % newton step

dx(:,1:4) = dx(:,1:4).*[n1>0, n2>0, n1s>0, n2s>0];

   x = x + dx;     
  
%fprintf('%4i  %6.2e \n', [it,  norm(dx)]);   %report difference

err = norm(dx);

end

% Go back to function iteration if newton does not converge

mu1H = 1./x(:,1); 
mu2H = 1./x(:,2); 
mu1F = 1./x(:,3); 
mu2F = 1./x(:,4);

Pi   = x(:,5).^(1/(gamma-1));

xold = [mu1H, mu2H, mu1F, mu2F, Pi];     
  
    it = 0; 
    
while err > 1e-7 && it < 500
   
    it = it + 1; 
    
mu1H = xold(:,1); 
mu2H = xold(:,2); 
mu1F = xold(:,3); 
mu2F = xold(:,4); 
Pi   = xold(:,5); 
    
    mu1H = (1 - 1/gamma + (mu1H./Pi).^(1-gamma).*(1/gamma-1/theta).*ui^(gamma-1)).^(-1).*(n1>0) + (n1==0);   
    mu2H = (max(1 - 1/gamma + (mu2H./Pi).^(1-gamma).*(1/gamma-1/theta).*(ui*ebar)^(gamma-1).*(n2.*xi + (1-xi)), 1-1/theta)).^(-1).*(n2>0) + (n2==0);
    mu1F = (1 - 1/gamma + (mu1F./Pi).^(1-gamma).*(1/gamma-1/theta).*(1/tau)^(gamma-1)).^(-1).*(n1s>0) + (n1s==0);   
    mu2F = (max(1 - 1/gamma + (mu2F./Pi).^(1-gamma).*(1/gamma-1/theta).*(ebar/tau)^(gamma-1).*(n2s.*xis + (1-xis)), 1-1/theta)).^(-1).*(n2s>0) + (n2s==0);
    Pi   = (n1.*(mu1H/ui).^(1-gamma) + n2.*(mu2H/ui/ebar).^(1-gamma) + ...
            n1s.*(tau*mu1F).^(1-gamma) + n2s.*(tau*mu2F/ebar).^(1-gamma)).^(1/(1-gamma));

xnew =  [mu1H, mu2H, mu1F, mu2F, Pi];     

%fprintf('%4i %6.2e \n',[it, norm(xnew - xold)]);   %report difference

err = norm(xnew-xold); 

xold = 1/20*xnew + 19/20*xold; 

end

mu1Hall(:,itt) = mu1H;
mu2Hall(:,itt) = mu2H; 
mu1Fall(:,itt) = mu1F;
mu2Fall(:,itt) = mu2F; 
Piall(:,itt)   = Pi; 

end
    
% Compute aggregate statistics

mu1H = mu1Hall; 
mu2H = mu2Hall; 
mu1F = mu1Fall;
mu2F = mu2Fall; 
Pi   = Piall; 

n1s  = [n1all(S/2+1:S,:); n1all(1:S/2,:)];
n2s  = [n2all(S/2+1:S,:); n2all(1:S/2,:)];

n1   = n1all;
n2   = n2all; 

a1H  = repmat(typei(:,3)', S, 1); 
a2H  = repmat(typei(:,3)', S, 1)*ebar; 
a1F  = ones(S, K); 
a2F  = ones(S, K)*ebar; 


p1H = mu1H./a1H; 
p2H = mu2H./a2H;
p1F = mu1F./a1F;
p2F = mu2F./a2F;


% Vectorize: 

n1   = n1(:); 
n2   = n2(:); 
n1s  = n1s(:); 
n2s  = n2s(:); 

mu1H = mu1H(:); 
mu2H = mu2H(:);
mu1F = mu1F(:);
mu2F = mu2F(:);

p1H  = p1H(:); 
p2H  = p2H(:); 
p1F  = p1F(:); 
p2F  = p2F(:); 

a1H  = a1H(:); 
a2H  = a2H(:);
a1F  = a1F(:);
a2F  = a2F(:); 

Pi   = Pi(:); 

wi   = kron(probi, ones(S,1)); 
wi   = wi/sum(wi); 

P   = (wi'*Pi.^(1-theta))^(1/(1-theta)); 

mui = (n1./mu1H.*(p1H./Pi).^(1-gamma) + n2./mu2H.*(p2H./Pi).^(1-gamma)+...
      n1s./mu1F.*(tau.*p1F./Pi).^(1-gamma) + n2s./mu2F.*(tau.*p2F./Pi).^(1-gamma)).^(-1);

Ai  = (n1.*(mu1H./mui).^(-gamma).*(a1H).^(gamma-1) + n2.*(mu2H./mui).^(-gamma).*a2H.^(gamma-1) + ...
       n1s.*(mu1F./mui).^(-gamma).*(a1F./tau).^(gamma-1) + n2s.*(mu2F./mui).^(-gamma).*(a2F./tau).^(gamma-1)).^(1/(gamma-1)); 

mu  = (wi'*(1./mui.*(Pi./P).^(1-theta))).^(-1);
   
A   = (wi'*((mui./mu).^(-theta).*Ai.^(theta-1))).^(1/(theta-1));
Y   = A;     


y1H = (p1H./Pi).^(-gamma).*(Pi/P).^(-theta).*Y; 
y2H = (p2H./Pi).^(-gamma).*(Pi/P).^(-theta).*Y; 
y1F = (tau*p1F./Pi).^(-gamma).*(Pi/P).^(-theta).*Y; 
y2F = (tau*p2F./Pi).^(-gamma).*(Pi/P).^(-theta).*Y; 


junk = reshape(y1F, S, K); 
y1Hs = [junk(S/2+1:S,:); junk(1:S/2,:)]; % uses symmetry in n draws
y1Hs = y1Hs(:); 

junk = reshape(y2F, S, K); 
y2Hs = [junk(S/2+1:S,:); junk(1:S/2,:)]; % uses symmetry in n draws
y2Hs = y2Hs(:); 

junk = reshape(p1F, S, K); 
p1Hs = [junk(S/2+1:S,:); junk(1:S/2,:)]; % uses symmetry in n draws
p1Hs = p1Hs(:); 

junk = reshape(p2F, S, K); 
p2Hs = [junk(S/2+1:S,:); junk(1:S/2,:)]; % uses symmetry in n draws
p2Hs = p2Hs(:); 

L   = wi'*(n1.*y1H./a1H + n2.*y2H./a2H + n1s.*tau.*y1F./a1F + n2s.*tau.*y2F./a2F);

% Naive Armington Elasticity (fix measures of producers)

lambdai    = n1.*(p1H./Pi).^(1-gamma) + n2.*(p2H./Pi).^(1-gamma);
si         = (Pi./P).^(1-theta);
lambda     = wi'*(si.*lambdai);
aweight    = wi'*(si.*lambdai/lambda.*(1-lambdai)/(1-lambda));
naiveArm   = gamma*aweight + theta*(1-aweight);

if abs(lambda-1) < 1e-3   % don't divide by 0
  
    aweight = 1;
    naiveArm = gamma;
    
end

impshare = 1 - lambdai; 
expshare = (n1.*p1Hs.*y1Hs.*tau + n2.*p2Hs.*y2Hs.*tau)./(n1.*p1H.*y1H + n2.*p2H.*y2H + n1.*p1Hs.*y1Hs.*tau + n2.*p2Hs.*y2Hs.*tau);


lambdai  = lambdai(:); 
impshare = impshare(:); 
expshare = expshare(:); 

ttrade   = impshare + expshare;

GLi      = 1 - abs(impshare - expshare)./ttrade;  

GL       = wi'*(si.*GLi);          % 1 means all trade is within industry

impshare = wi'*(si.*impshare); 

% Correct Armington Elasticity (trade costs affect markups and measures of
% producers

if comp_armington

params = [gamma; theta; ebar; S; K];
    
step   = 1.005; 
taunew = tau*step;

impsharenew = armington(taunew, mu1Hall, mu2Hall, mu1Fall, mu2Fall, n1all, n2all, typei, probi, params); 

Arm = log(impsharenew./(1-impsharenew)/(impshare/(1-impshare)))/log(step);

else
    
    Arm = 0;
    
end

% losses from misallocation

Aibest  = (n1.*a1H.^(gamma-1) + n2.*a2H.^(gamma-1) + ...
           n1s.*(a1F./tau).^(gamma-1) + n2s.*(a2F/tau).^(gamma-1)).^(1/(gamma-1)); 
   
Abest   = (wi'*(Aibest.^(theta-1))).^(1/(theta-1));

meanmu = wi'*(n1.*mu1H + n2.*mu2H + n1s.*mu1F + n2s.*mu2F)./(wi'*(n1 + n2 + n1s + n2s)); 

fprintf('\n');
display('Equilibrium Implications');
fprintf('\n');
fprintf('Y                   = %7.3f \n',  Y);  
fprintf('A                   = %7.3f \n',  A);
fprintf('L                   = %7.3f \n',  L);
fprintf('aggregate profits   = %7.3f \n',  (P*Y - L)/P); 
fprintf('A loss, *100        = %7.3f \n',  log(Abest/A)*100);  
fprintf('mean markup         = %7.3f \n',  meanmu);
fprintf('aggregate markup    = %7.3f \n',  P*A);
fprintf('import share        = %7.3f \n',  impshare);
fprintf('trade elastic       = %7.3f \n',  -Arm);

fprintf('\n');
