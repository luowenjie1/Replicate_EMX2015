function [y, yd] = solveequil(x, tau, ui, xi, xis, n1, n2, n1s, n2s, gamma, theta, ebar)

S = size(x,1); 

y = zeros(S, 5); 

a = zeros(4,1); 

a(1) = ui; 
a(2) = ui*ebar;
a(3) = 1/tau; 
a(4) = ebar/tau; 

y(:,1) = x(:,1) - (1-1/gamma) - x(:,1).^(gamma-1).*x(:,5)*(1/gamma-1/theta).*a(1).^(gamma-1);
y(:,2) = x(:,2) - (1-1/gamma) - x(:,2).^(gamma-1).*x(:,5)*(1/gamma-1/theta).*a(2).^(gamma-1).*(n2.*xi + (1-xi));
y(:,3) = x(:,3) - (1-1/gamma) - x(:,3).^(gamma-1).*x(:,5)*(1/gamma-1/theta).*a(3).^(gamma-1);
y(:,4) = x(:,4) - (1-1/gamma) - x(:,4).^(gamma-1).*x(:,5)*(1/gamma-1/theta).*a(4).^(gamma-1).*(n2s.*xis + (1-xis));
y(:,5) = n1.*(x(:,1).*a(1)).^(gamma-1) + n2.*(x(:,2).*a(2)).^(gamma-1) + ...
         n1s.*(x(:,3).*a(3)).^(gamma-1) + n2s.*(x(:,4).*a(4)).^(gamma-1) - 1./x(:,5); 

y(:,1) = y(:,1).*(n1>0); % don't bother if n1=0 since it doesn't matter
y(:,2) = y(:,2).*(n2>0); % don't bother if n1=0 since it doesn't matter
y(:,3) = y(:,3).*(n1s>0); % don't bother if n1=0 since it doesn't matter
y(:,4) = y(:,4).*(n2s>0); % don't bother if n1=0 since it doesn't matter
         
yd = zeros(S, 5, 5); 
     
yd(:,1,1) = 1 - (gamma-1)*x(:,1).^(gamma-2).*x(:,5)*(1/gamma-1/theta).*a(1).^(gamma-1); 
yd(:,2,2) = 1 - (gamma-1)*x(:,2).^(gamma-2).*x(:,5)*(1/gamma-1/theta).*a(2).^(gamma-1).*(n2.*xi + (1-xi));
yd(:,3,3) = 1 - (gamma-1)*x(:,3).^(gamma-2).*x(:,5)*(1/gamma-1/theta).*a(3).^(gamma-1);
yd(:,4,4) = 1 - (gamma-1)*x(:,4).^(gamma-2).*x(:,5)*(1/gamma-1/theta).*a(4).^(gamma-1).*(n2s.*xis + (1-xis));

yd(:,1,5) = x(:,1).^(gamma-1).*(1/gamma-1/theta).*a(1).^(gamma-1);
yd(:,2,5) = x(:,2).^(gamma-1).*(1/gamma-1/theta).*a(2).^(gamma-1).*(n2.*xi + (1-xi));
yd(:,3,5) = x(:,3).^(gamma-1).*(1/gamma-1/theta).*a(3).^(gamma-1);
yd(:,4,5) = x(:,4).^(gamma-1).*(1/gamma-1/theta).*a(4).^(gamma-1).*(n2s.*xis + (1-xis));

yd(:,5,1) = (gamma-1).*n1.*(x(:,1).*a(1)).^(gamma-1)./x(:,1);
yd(:,5,2) = (gamma-1).*n2.*(x(:,2).*a(2)).^(gamma-1)./x(:,2); 
yd(:,5,3) = (gamma-1).*n1s.*(x(:,3).*a(3)).^(gamma-1)./x(:,3); 
yd(:,5,4) = (gamma-1).*n2s.*(x(:,4).*a(4)).^(gamma-1)./x(:,4); 
yd(:,5,5) = x(:,5).^(-2);

