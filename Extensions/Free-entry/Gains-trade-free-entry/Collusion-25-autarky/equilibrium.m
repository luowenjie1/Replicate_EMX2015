rand('state', 0) 

if 1

xi  = [1; 0];  
    
typei = gridmake(om1, om2, uu, xi); 
probi = gridmake(pom, wu, [lambda; 1-lambda]);
probi = probi(:,1).*probi(:,2).*probi(:,3); 

else

    xi = 0; 
    
typei = gridmake(om1, om2, uu, xi); 
probi = gridmake(pom, wu, 1);
probi = probi(:,1).*probi(:,2).*probi(:,3); 

end    

K = size(typei,1); 

n1all  = zeros(S, K); 
n2all  = n1all;

parfor  itt = 1 : K
    
n            = binornd(N, typei(itt,1), S, 1) + 1;  % there is at least one firm in an industry
n2all(:,itt) = binornd(n, typei(itt,2));
n1all(:,itt) = n - n2all(:,itt); 

end

clear n

mu1Hall = zeros(S, K);
mu2Hall = zeros(S, K); 
Piall   = zeros(S, K); 

fprintf('\n')
fprintf('Solve each industrys equilibrium \n')
fprintf('\n')

parfor itt = 1 : K
        
n1  = n1all(:,itt); 
n2  = n2all(:,itt); 
    
ui  = typei(itt,3); 
xi  = typei(itt,4);   % can only collude if n2>1
    
xi   = xi*(n2>1);       % collude dummy

% compute industry equilibrium for a given u(i), xi(i), xis(i)

mu1H = gamma/(gamma-1)*ones(S,1).*(n1>0) + (n1==0); 
mu2H = gamma/(gamma-1)*ones(S,1).*(n2>0) + (n2==0); 

Pi   = (n1.*(mu1H/ui).^(1-gamma)   + n2.*(mu2H/ui/ebar).^(1-gamma)).^(1/(1-gamma));

xold = [mu1H, mu2H, Pi];     

err = 1000; 

% A few function iterations

for it = 1: 25
   
mu1H = xold(:,1); 
mu2H = xold(:,2); 
Pi   = xold(:,3); 

    mu1H = (1 - 1/gamma + (mu1H./Pi).^(1-gamma).*(1/gamma-1/theta).*ui^(gamma-1)).^(-1).*(n1>0) + (n1==0);   
    mu2H = (max(1 - 1/gamma + (mu2H./Pi).^(1-gamma).*(1/gamma-1/theta).*(ui*ebar)^(gamma-1).*(n2.*xi + (1-xi)), 1-1/theta)).^(-1).*(n2>0) + (n2==0);
    Pi   = (n1.*(mu1H/ui).^(1-gamma) + n2.*(mu2H/ui/ebar).^(1-gamma)).^(1/(1-gamma));

xnew =  [mu1H, mu2H, Pi];     

%fprintf('%4i %6.2e \n',[it, norm(xnew - xold)]);   %report difference

err  = norm(xnew - xold);

xold = 1/20*xnew + 19/20*xold; 

end
   
x = [1./mu1H, 1./mu2H, Pi.^(gamma-1)];

% Newton

it = 0; 

while err > 1e-7 && it < 250
    
    it = it + 1; 
    
[y, yd] = solveequil(x, ui, xi,  n1, n2, gamma, theta, ebar);

dx = - arrayinv(y, yd); % newton step

dx(:,1:2) = dx(:,1:2).*[n1>0, n2>0];

   x = x + dx;     
  
%fprintf('%4i  %6.2e \n', [it,  norm(dx)]);   %report difference

err = norm(dx);

end

% Go back to function iteration if newton does not converge

mu1H = 1./x(:,1); 
mu2H = 1./x(:,2); 

Pi   = x(:,3).^(1/(gamma-1));

xold = [mu1H, mu2H, Pi];     
  
    it = 0; 
    
while err > 1e-7 && it < 500
   
    it = it + 1; 
    
mu1H = xold(:,1); 
mu2H = xold(:,2); 

Pi   = xold(:,3); 
    
    mu1H = (1 - 1/gamma + (mu1H./Pi).^(1-gamma).*(1/gamma-1/theta).*ui^(gamma-1)).^(-1).*(n1>0) + (n1==0);   
    mu2H = (max(1 - 1/gamma + (mu2H./Pi).^(1-gamma).*(1/gamma-1/theta).*(ui*ebar)^(gamma-1).*(n2.*xi + (1-xi)), 1-1/theta)).^(-1).*(n2>0) + (n2==0);
    Pi   = (n1.*(mu1H/ui).^(1-gamma) + n2.*(mu2H/ui/ebar).^(1-gamma)).^(1/(1-gamma));

xnew =  [mu1H, mu2H, Pi];     

%fprintf('%4i %6.2e \n',[it, norm(xnew - xold)]);   %report difference

err = norm(xnew-xold); 

xold = 1/20*xnew + 19/20*xold; 

end

mu1H(n1==0)  = gamma/(gamma-1); 
mu2H(n2==0)  = gamma/(gamma-1); 

mu1Hall(:,itt) = mu1H;
mu2Hall(:,itt) = mu2H; 
Piall(:,itt)   = Pi; 

end

% Compute aggregate Price level

mu1H = mu1Hall; 
mu2H = mu2Hall; 

Pi   = Piall; 

n1   = n1all;
n2   = n2all; 

a1H  = repmat(typei(:,3)', S, 1); 
a2H  = repmat(typei(:,3)', S, 1)*ebar; 

p1H = mu1H./a1H; 
p2H = mu2H./a2H;

% Vectorize: 

n1   = n1(:); 
n2   = n2(:);  

mu1H = mu1H(:); 
mu2H = mu2H(:);

p1H  = p1H(:); 
p2H  = p2H(:); 

a1H  = a1H(:); 
a2H  = a2H(:); 

Pi   = Pi(:); 

wi   = kron(probi, ones(S,1)); 
wi   = wi/sum(wi); 

P    = (wi'*Pi.^(1-theta))^(1/(1-theta)); 

mui  = (n1./mu1H.*(p1H./Pi).^(1-gamma) + n2./mu2H.*(p2H./Pi).^(1-gamma)).^(-1);

Ai   = (n1.*(mu1H./mui).^(-gamma).*(a1H).^(gamma-1) + n2.*(mu2H./mui).^(-gamma).*a2H.^(gamma-1)).^(1/(gamma-1)); 

mu   = (wi'*(1./mui.*(Pi./P).^(1-theta))).^(-1);
   
A    = (wi'*((mui./mu).^(-theta).*Ai.^(theta-1))).^(1/(theta-1));
Y    = A;    

Aibest  = (n1.*a1H.^(gamma-1) + n2.*a2H.^(gamma-1)).^(1/(gamma-1)); 
   
Abest   = (wi'*(Aibest.^(theta-1))).^(1/(theta-1));


a1H  = repmat(typei(:,3)', S, 1); 
a2H  = repmat(typei(:,3)', S, 1)*ebar; 


% Next, equilibrium if n1 = n1 + 1 (profits from entering low-type in
% Home): call this c1 (case 1)

mu1Hc1 = zeros(S, K);
mu2Hc1 = zeros(S, K); 
Pic1   = zeros(S, K); 

fprintf('\n')
fprintf('Solve each industrys equilibrium assuming n1 = n1 + 1 \n')
fprintf('\n')

parfor itt = 1 : K
        
n1  = n1all(:,itt); 
n2  = n2all(:,itt); 
    
ui  = typei(itt,3); 
xi  = typei(itt,4);   % can only collude if n2>1
          
n1  = n1 + 1; 

xi   = xi*(n2>1);       % collude dummy

% use original equilibrium as guess

mu1H = mu1Hall(:,itt);
mu2H = mu2Hall(:,itt);
 
%Pi   = Piall(:,itt); 
Pi   = (n1.*(mu1H/ui).^(1-gamma) + n2.*(mu2H/ui/ebar).^(1-gamma)).^(1/(1-gamma));

xold = [mu1H, mu2H,  Pi];     

err = 1000; 

% A few function iterations

for it = 1: 25
   
mu1H = xold(:,1); 
mu2H = xold(:,2); 

Pi   = xold(:,3); 

    mu1H = (1 - 1/gamma + (mu1H./Pi).^(1-gamma).*(1/gamma-1/theta).*ui^(gamma-1)).^(-1).*(n1>0) + (n1==0);   
    mu2H = (max(1 - 1/gamma + (mu2H./Pi).^(1-gamma).*(1/gamma-1/theta).*(ui*ebar)^(gamma-1).*(n2.*xi + (1-xi)), 1-1/theta)).^(-1).*(n2>0) + (n2==0);
    Pi   = (n1.*(mu1H/ui).^(1-gamma) + n2.*(mu2H/ui/ebar).^(1-gamma)).^(1/(1-gamma));

xnew =  [mu1H, mu2H, Pi];     

%fprintf('%4i %6.2e \n',[it, norm(xnew - xold)]);   %report difference

err  = norm(xnew - xold);

xold = 1/20*xnew + 19/20*xold; 

end
   
x = [1./mu1H, 1./mu2H, Pi.^(gamma-1)];

% Newton

it = 0; 

while err > 1e-7 && it < 250
    
    it = it + 1; 
    
[y, yd] = solveequil(x, ui, xi, n1, n2, gamma, theta, ebar);

dx = - arrayinv(y, yd); % newton step

dx(:,1:2) = dx(:,1:2).*[n1>0, n2>0];

   x = x + dx;     
  
%fprintf('%4i  %6.2e \n', [it,  norm(dx)]);   %report difference

err = norm(dx);

end

% Go back to function iteration if newton does not converge

mu1H = 1./x(:,1); 
mu2H = 1./x(:,2); 

Pi   = x(:,3).^(1/(gamma-1));

xold = [mu1H, mu2H, Pi];     
  
    it = 0; 
    
while err > 1e-7 && it < 500
   
    it = it + 1; 
    
mu1H = xold(:,1); 
mu2H = xold(:,2); 
Pi   = xold(:,3); 
    
    mu1H = (1 - 1/gamma + (mu1H./Pi).^(1-gamma).*(1/gamma-1/theta).*ui^(gamma-1)).^(-1).*(n1>0) + (n1==0);   
    mu2H = (max(1 - 1/gamma + (mu2H./Pi).^(1-gamma).*(1/gamma-1/theta).*(ui*ebar)^(gamma-1).*(n2.*xi + (1-xi)), 1-1/theta)).^(-1).*(n2>0) + (n2==0);
    Pi   = (n1.*(mu1H/ui).^(1-gamma) + n2.*(mu2H/ui/ebar).^(1-gamma)).^(1/(1-gamma));

xnew =  [mu1H, mu2H,  Pi];     

%fprintf('%4i %6.2e \n',[it, norm(xnew - xold)]);   %report difference

err = norm(xnew-xold); 

xold = 1/20*xnew + 19/20*xold; 

end

mu1Hc1(:,itt) = mu1H;
mu2Hc1(:,itt) = mu2H; 
Pic1(:,itt)   = Pi; 

end

a1H  = repmat(typei(:,3)', S, 1);  

profit1H = (mu1Hc1 - 1)./a1H.*(mu1Hc1./a1H./Pic1).^(-gamma).*(Pic1./P).^(-theta)*Y;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Next, equilibrium if n2 = n2 + 1 (profits from entering high-type in
% Home): call this c2 (case 2)

mu1Hc2 = zeros(S, K);
mu2Hc2 = zeros(S, K); 
Pic2   = zeros(S, K); 

fprintf('\n')
fprintf('Solve each industrys equilibrium assuming n2 = n2 + 1 \n')
fprintf('\n')

for itt = 1 : K

n1  = n1all(:,itt); 
n2  = n2all(:,itt); 
    
ui  = typei(itt,3); 
xi  = typei(itt,4);   % can only collude if n2>1
          
n2  = n2 + 1; 

xi   = xi*(n2>1);       % collude dummy

% use original equilibrium as guess

mu1H = mu1Hall(:,itt);
mu2H = mu2Hall(:,itt); 
%Pi   = Piall(:,itt); 
Pi   = (n1.*(mu1H/ui).^(1-gamma) + n2.*(mu2H/ui/ebar).^(1-gamma)).^(1/(1-gamma));

xold = [mu1H, mu2H, Pi];     

err = 1000; 

% A few function iterations

for it = 1: 25
   
mu1H = xold(:,1); 
mu2H = xold(:,2); 
Pi   = xold(:,3); 

    mu1H = (1 - 1/gamma + (mu1H./Pi).^(1-gamma).*(1/gamma-1/theta).*ui^(gamma-1)).^(-1).*(n1>0) + (n1==0);   
    mu2H = (max(1 - 1/gamma + (mu2H./Pi).^(1-gamma).*(1/gamma-1/theta).*(ui*ebar)^(gamma-1).*(n2.*xi + (1-xi)), 1-1/theta)).^(-1).*(n2>0) + (n2==0);
    Pi   = (n1.*(mu1H/ui).^(1-gamma) + n2.*(mu2H/ui/ebar).^(1-gamma)).^(1/(1-gamma));

xnew =  [mu1H, mu2H, Pi];     

%fprintf('%4i %6.2e \n',[it, norm(xnew - xold)]);   %report difference

err  = norm(xnew - xold);

xold = 1/20*xnew + 19/20*xold; 

end
   
x = [1./mu1H, 1./mu2H, Pi.^(gamma-1)];

% Newton

it = 0; 

while err > 1e-7 && it < 250
    
    it = it + 1; 
    
[y, yd] = solveequil(x, ui, xi, n1, n2, gamma, theta, ebar);

dx = - arrayinv(y, yd); % newton step

dx(:,1:2) = dx(:,1:2).*[n1>0, n2>0];

   x = x + dx;     
  
%fprintf('%4i  %6.2e \n', [it,  norm(dx)]);   %report difference

err = norm(dx);

end

% Go back to function iteration if newton does not converge

mu1H = 1./x(:,1); 
mu2H = 1./x(:,2); 

Pi   = x(:,3).^(1/(gamma-1));

xold = [mu1H, mu2H, Pi];     
  
    it = 0; 
    
while err > 1e-7 && it < 500
   
    it = it + 1; 
    
mu1H = xold(:,1); 
mu2H = xold(:,2); 
Pi   = xold(:,3); 
    
    mu1H = (1 - 1/gamma + (mu1H./Pi).^(1-gamma).*(1/gamma-1/theta).*ui^(gamma-1)).^(-1).*(n1>0) + (n1==0);   
    mu2H = (max(1 - 1/gamma + (mu2H./Pi).^(1-gamma).*(1/gamma-1/theta).*(ui*ebar)^(gamma-1).*(n2.*xi + (1-xi)), 1-1/theta)).^(-1).*(n2>0) + (n2==0);
    Pi   = (n1.*(mu1H/ui).^(1-gamma) + n2.*(mu2H/ui/ebar).^(1-gamma)).^(1/(1-gamma));

xnew =  [mu1H, mu2H, Pi];     

%fprintf('%4i %6.2e \n',[it, norm(xnew - xold)]);   %report difference

err = norm(xnew-xold); 

xold = 1/20*xnew + 19/20*xold; 

end

mu1Hc2(:,itt) = mu1H;
mu2Hc2(:,itt) = mu2H; 
Pic2(:,itt)   = Pi; 

end

a2H  = repmat(typei(:,3)', S, 1)*ebar; 

profit2H = (mu2Hc2 - 1)./a2H.*(mu2Hc2./a2H./Pic2).^(-gamma).*(Pic2./P).^(-theta)*Y;


% Next, compute expected profits
% Vectorize: 

profit1H  = profit1H(:); % S*K
profit2H  = profit2H(:); % S*K

wi   = kron(probi, ones(S,1));   % probability of entering any one of the S*K sectors
wi   = wi/sum(wi); 

omega1 = repmat(typei(:,1)', S, 1); 
omega2 = repmat(typei(:,2)', S, 1); 
omega1 = omega1(:); 
omega2 = omega2(:); 

n1   = n1all(:);
n2   = n2all(:); 

Eprofits = wi'*(omega1.*(1-omega2).*(profit1H) + ...
                omega1.*omega2.*(profit2H))/P; 
            
Aprofits = wi'*(n1.*(profit1H) + n2.*(profit2H))/P/N;
            
fprintf('\n');
fprintf('Average real profits (of an incumbent) (*100)  = %7.3f \n',  (P*Y - L)/P/N*100);
fprintf('Average real profits (of an entrant) (*100)    = %7.3f \n',  Aprofits*100);
fprintf('Expected real profits (*100)                   = %7.3f \n',  Eprofits*100);

load fixedcost

fprintf('\n');
fprintf('\n');
display('Equilibrium Implications');
fprintf('\n');
fprintf('Y                            = %7.3f \n',  Y);  
fprintf('A                            = %7.3f \n',  A);
fprintf('L                            = %7.3f \n',  L);
fprintf('A loss, *100                 = %7.3f \n',  log(Abest/A)*100);  
fprintf('Total fixed costs            = %7.3f \n',  fixedcost*N);
fprintf('C                            = %7.3f \n',  Y - fixedcost*N);
fprintf('aggregate profits            = %7.3f \n',  (P*Y - L)/P); 
fprintf('aggregate markup             = %7.3f \n',  P*A);
fprintf('profits - fixed cost         = %7.3f \n',  (P*Y - L)/P - fixedcost*N); 


fprintf('E profits vs. fixed cost     = %7.3f %7.3f \n', [Eprofits*100, fixedcost*100]); 

