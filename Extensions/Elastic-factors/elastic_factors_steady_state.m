function [C,K,Y,L,U] = elastic_factors_steady_state(parameters,A,M);

% steady state of model with capital accumulation, elastic labor supply

alpha = parameters(1) ; %% output elasticity of capital in aggregate production function
delta = parameters(2) ; %% capital depreciation
beta  = parameters(3) ; %% discount factor
chi   = parameters(4) ; %% weight on leisure (chosen so L=1)
eta   = parameters(5) ; %% inverse Frisch elasticity labor supply

rho   = 1/beta - 1;

KtoL = ((alpha/(rho+delta))*(A/M))^(1/(1-alpha)); %% K/L
YtoL = A*(KtoL)^alpha;                            %% Y/L
KtoY = KtoL/YtoL;                                 %% K/Y
CtoY = 1 - delta*KtoY;                            %% C/Y

L    = ((1/chi)*(1-alpha)/(M-(alpha*delta)/(rho+delta)))^(1/(eta+1)); %% L

K    = KtoL*L;
Y    = YtoL*L;
C    = CtoY*Y;


U    = log(C) - chi*(L^(1+eta))/(1+eta);



