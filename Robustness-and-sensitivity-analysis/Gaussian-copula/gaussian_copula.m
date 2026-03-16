function [ux,uy] = gaussian_copula(N,taurho)
%GAUSSIAN_COPULA ux,uy are draws of length N from a gaussian copula with
%"Kendall's tau" = taurho

uu = nodeunif(N,eps^(1/3),1-eps^(1/3));
ux = uu(randperm(N));
uy = uu(randperm(N));

theta = sin(taurho*pi/2); % linear correlation coefficient corresponding to taurho   
    
u1 = norminv(ux);
u2 = norminv(uy);
    
z1 = u1;
z2 = u1*theta+u2*sqrt(1-theta^2);

ux = normcdf(z1);
uy = normcdf(z2);

end

