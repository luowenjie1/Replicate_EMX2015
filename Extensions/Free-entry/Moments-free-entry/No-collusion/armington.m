function impshare = armington(tau, mu1Hall, mu2Hall, mu1Fall, mu2Fall, n1all, n2all, typei, probi, params)

gamma = params(1);
theta = params(2);
ebar  = params(3);
S     = params(4);
K     = params(5);


parfor itt = 1 : K
        
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

mu1H = mu1Hall(:,itt); 
mu2H = mu2Hall(:,itt); 
mu1F = mu1Fall(:,itt); 
mu2F = mu2Fall(:,itt);

Pi   = (n1.*(mu1H/ui).^(1-gamma)   + n2.*(mu2H/ui/ebar).^(1-gamma) + ...
        n1s.*(tau*mu1F).^(1-gamma) + n2s.*(tau*mu2F/ebar)).^(1/(1-gamma));

xold = [mu1H, mu2H, mu1F, mu2F, Pi];     

err = 1000; 

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

err = norm(xnew - xold);

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
Pi   = Piall; 

n1   = n1all;
n2   = n2all; 

clear n1all n2all mu1Hall mu2Hall mu1Fall mu2Fall;

a1H  = repmat(typei(:,3)', S, 1); 
a2H  = repmat(typei(:,3)', S, 1)*ebar; 


p1H = mu1H./a1H; 
p2H = mu2H./a2H;

% Vectorize: 

n1   = n1(:); 
n2   = n2(:); 

p1H  = p1H(:); 
p2H  = p2H(:); 

Pi   = Pi(:); 

wi   = kron(probi, ones(S,1)); 
wi   = wi/sum(wi); 

P   = (wi'*Pi.^(1-theta))^(1/(1-theta)); 

% Naive Armington Elasticity (fix measures of producers)

lambdai    = n1.*(p1H./Pi).^(1-gamma) + n2.*(p2H./Pi).^(1-gamma);
si         = (Pi./P).^(1-theta);
lambda     = wi'*(si.*lambdai);

impshare   = 1 - lambda; 