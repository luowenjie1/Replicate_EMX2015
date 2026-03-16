function [y, yd] = solveequil(x, ui, xi, n1, n2, gamma, theta, ebar)

S = size(x,1); 

y = zeros(S, 3); 

a = zeros(2,1); 

a(1) = ui; 
a(2) = ui*ebar;


y(:,1) = x(:,1) - (1-1/gamma) - x(:,1).^(gamma-1).*x(:,3)*(1/gamma-1/theta).*a(1).^(gamma-1);
y(:,2) = x(:,2) - (1-1/gamma) - x(:,2).^(gamma-1).*x(:,3)*(1/gamma-1/theta).*a(2).^(gamma-1).*(n2.*xi + (1-xi));
y(:,3) = n1.*(x(:,1).*a(1)).^(gamma-1) + n2.*(x(:,2).*a(2)).^(gamma-1) - 1./x(:,3); 

y(:,1) = y(:,1).*(n1>0); % don't bother if n1=0 since it doesn't matter
y(:,2) = y(:,2).*(n2>0); % don't bother if n1=0 since it doesn't matter
         
yd = zeros(S, 3, 3); 
     
yd(:,1,1) = 1 - (gamma-1)*x(:,1).^(gamma-2).*x(:,3)*(1/gamma-1/theta).*a(1).^(gamma-1); 
yd(:,2,2) = 1 - (gamma-1)*x(:,2).^(gamma-2).*x(:,3)*(1/gamma-1/theta).*a(2).^(gamma-1).*(n2.*xi + (1-xi));

yd(:,1,3) = x(:,1).^(gamma-1).*(1/gamma-1/theta).*a(1).^(gamma-1);
yd(:,2,3) = x(:,2).^(gamma-1).*(1/gamma-1/theta).*a(2).^(gamma-1).*(n2.*xi + (1-xi));

yd(:,3,1) = (gamma-1).*n1.*(x(:,1).*a(1)).^(gamma-1)./x(:,1);
yd(:,3,2) = (gamma-1).*n2.*(x(:,2).*a(2)).^(gamma-1)./x(:,2); 
yd(:,3,3) = x(:,3).^(-2);

