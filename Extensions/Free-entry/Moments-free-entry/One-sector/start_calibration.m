clear all
clc

%%%%% parameters to calibrate

N       = 176;                  % number of firms that draw
om1     = 0.0972513900; 
om2     = 0.1734328984; 
ebar    = 1.838675631;          % productivity of type 2
su      = 0.1469919780;         % standard deviation of Home/Foreign productivity gap
tau     = 1.1339970969; 



x      = [N, om1, om2, ebar, su, tau];
nvars  = length(x);

%%%%% constraints on x

lb     = [100,   0.01,  0.01,     1.1,  0.05,  1.01];        % lower bounds
ub     = [200,   0.2,   0.2 ,     3     ,  0.7  ,  2    ];        % upper bounds

intcon = 1;   %  x(1) must be an integer, number of firms

warning('off', 'MATLAB:nearlySingularMatrix')

%%%%%% settings for genetic algorithm

A   = []; b   = [];
Aeq = []; beq = []; nonlcon = [];

addpath(genpath('/scratch/jab772/compecom')); 
matlabpool open 12

%gaoptions = gaoptimset(gaoptions,'UseParallel','always','Display','iter','InitialPopulation',xstart);
gaoptions = gaoptimset('UseParallel', 'always', 'Display','iter','InitialPopulation',x);
x = ga(@objective, nvars, A, b, Aeq, beq, lb, ub, [], intcon, gaoptions);

%options = gaoptimset('Generations',100,'Vectorized','on','UseParallel','always');    
%[x,fval,exitflag,output] = ga(@(x)objective(x,params),nvars,A,b,Aeq,beq,lb,ub,nonlcon,intcon);

xga=x;

save Results xga

matlabpool close

x     = x';
xsave = x;   % save in case fminsearch goes wild

x = fminsearch('objective', x, optimset('Display', 'iter')); 

save Results xsave xga

