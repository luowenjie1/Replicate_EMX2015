% -------------------------------------------------------------------------
% prod_dlwg.m
% Purpose: Estimate sector-level gross-output production coefficients
%          (DLW-style proxy approach) for sectors 10-31.
%
% [DATA AVAILABILITY NOTICE]
% Requires confidential MATLAB files in ../../../production/ named datag_10.mat
% ... datag_31.mat, generated from restricted plant-level micro data.
% -------------------------------------------------------------------------

clc
clear all
cd('../../../production')

global lnqh lagqh lnk lagk lnl lagl lnm lagm

% Store coefficients by 2-digit sector: [sic2, 9 parameters].
prod_est = zeros(22,10);
for j = 10:1:31
    a = ['datag_' num2str(j) '.mat'];
    load(a)

    % OLS starting values from second-order translog in labor/capital/material.
    xmat = [ones(size(lnqh,1),1) lnl lnk lnm lnl.^2 lnk.^2 lnm.^2 lnl.*lnk lnl.*lnm lnk.*lnm];
    est_OLS = (xmat'*xmat)\(xmat'*lnqh);

    % Nonlinear objective in obj_dlwg.m (GMM-style orthogonality condition).
    options = optimset('MaxFunEvals', 1e+10);
    est = fminsearch(@obj_dlwg, est_OLS(2:end)', options);

    disp('input elasticity')
    [est;est_OLS(2:end)']

    al  = est(1);
    ak  = est(2);
    am  = est(3);
    al2 = est(4);
    ak2 = est(5);
    am2 = est(6);
    alk = est(7);
    alm = est(8);
    akm = est(9);

    % Inspect implied elasticity distributions.
    prctile(am+2*am2*lnm+akm*lnk+alm*lnl, [5 10 50 90 95])
    prctile(al+2*al2*lnl+alm*lnm+alk*lnk, [5 10 50 90 95])

    prod_est(j-9,:) = [j est];
end

% Export coefficients for Stata scripts.
csvwrite('prod_est_g.csv', prod_est);
