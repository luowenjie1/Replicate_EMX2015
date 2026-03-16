function out = obj_dlwg(in)
% -------------------------------------------------------------------------
% obj_dlwg.m
% Objective used by prod_dlwg.m for production-function estimation.
% Minimizes the quadratic form of sample moments E[z_t * xi_t] where xi_t is
% the innovation in latent productivity omega_t after projecting on lagged
% productivity terms.
% -------------------------------------------------------------------------

global lnqh lagqh lnk lagk lnl lagl lnm lagm

al  = in(1);
ak  = in(2);
am  = in(3);
al2 = in(4);
ak2 = in(5);
am2 = in(6);
alk = in(7);
alm = in(8);
akm = in(9);

% Current and lagged productivity residuals implied by candidate parameters.
omega  = lnqh  - al*lnl  - al2*lnl.^2  - ak*lnk  - ak2*lnk.^2  - alk*lnl.*lnk  - am*lnm  - am2*lnm.^2  - alm*lnl.*lnm  - akm*lnk.*lnm;
lomega = lagqh - al*lagl - al2*lagl.^2 - ak*lagk - ak2*lagk.^2 - alk*lagl.*lagk - am*lagm - am2*lagm.^2 - alm*lagl.*lagm - akm*lagk.*lagm;

nobs = size(omega,1);

% Project omega on polynomial in lagged omega; xi is innovation.
xmat = [ones(nobs,1) lomega lomega.^2];
xi   = omega - xmat*((xmat'*xmat)\(xmat'*omega));

% Instrument set (lagged inputs and interactions).
z = [ones(nobs,1) lagl lnk lagm lagl.^2 lnk.^2 lagm.^2 lnk.*lagl lagl.*lagm lnk.*lagm];

% Quadratic criterion: || Z' * xi ||^2.
out = (z'*xi)'*(z'*xi);
end
