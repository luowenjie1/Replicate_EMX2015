clc
clear all
cd ('C:\Users\yx42\Documents\Research Projects\Competition and Markup\Revision\production')

global lnqh lagqh lnk lagk lnl lagl lnm lagm

prod_est=zeros(22,10);
for j=10:1:31
   a=['datag_' num2str(j) '.mat'];
    load(a)
%OLS prediction
xmat=[ones(size(lnqh,1),1) lnl  lnk lnm lnl.^2 lnk.^2 lnm.^2 lnl.*lnk lnl.*lnm lnk.*lnm ];
est_OLS=(xmat'*xmat)\(xmat'*lnqh);
%disp('OLS estimates')
%[est_OLS]

options = optimset('MaxFunEvals', 1e+10);  %'TolX',1e-10,'TolFun', 1e-6,'
est=fminsearch(@obj_dlwg, est_OLS(2:end)',options);
disp('input elasticity')
[est;est_OLS(2:end)']


al=est(1);
ak=est(2);
am=est(3);
al2=est(4);
ak2=est(5);
am2=est(6);
alk=est(7);
alm=est(8);
akm=est(9);

prctile(am+2*am2*lnm+akm*lnk+alm*lnl, [5 10 50 90 95])
prctile(al+2*al2*lnl+alm*lnm+alk*lnk, [5 10 50 90 95])

%save estimates
prod_est(j-9,:)=[j est];

end

csvwrite('prod_est_g.csv', prod_est);
