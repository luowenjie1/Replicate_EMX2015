function [ux,uy] = gumbel_copula(N,taurho)
%GUMBEL_COPULA ux,uy are draws of length N from a gumbel copula with
%"Kendall's tau" = taurho

uu = nodeunif(N,eps^(1/3),1-eps^(1/3));
ux = uu(randperm(N));
uy = uu(randperm(N));

% Gumbel's copula

sold  = ones(size(ux))*eps;

rho   = 1/(1-taurho);

for it = 1:200

s = sold;

Ys  = log(s).*s - rho*(s-ux);   % equation implicitly defines s given ux
Yps = 1 + log(s) - rho;         % derivative

sold = s - 1/2*Yps.^(-1).*Ys;

if norm(s-sold)<1e-10, break, end

end

if norm(Ys)>1e-5
       disp('Warning: Gumbel has not converged')
end

data = [exp(log(s).*uy.^(1/rho)),exp(log(s).*(1-uy).^(1/rho))];

ux = data(:,1); 
uy = data(:,2);


end

