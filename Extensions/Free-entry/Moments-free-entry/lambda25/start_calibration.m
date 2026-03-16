clear all
clc

%%%%% parameters to calibrate

N       = 160;                                  % number of firms that draw
ebar    = 1.7176;                          % productivity of type 2
su      = 0.1737;                         % standard deviation of Home/Foreign productivity gap
tau     = 1.1356; 

om1     = [0.1583 , 0.0116, 0.3131]; 
om2     = [0.1535, 0.0587,0.2596]; 

pom     = [0.4681, 0.3000, 0.4239, ...
           0.5701,0.6364 ,0.6787 , ...
           0.4236, 0.7156, 0.6289];

x      = [N, ebar, su, tau, om1, om2,pom];

nvars  = length(x);

%%%%% constraints on x

lb     = [75,   1.1,  0.05,  1.01,...
          0.01,  0.01,0.01,0.01,  0.01,0.01,...
          zeros(1,9)];        % lower bounds
ub     = [200,     3,  0.7 ,  2   ,...
          0.5,   0.5,   0.5,0.5,   0.5,   0.5,...
          ones(1,9)];        % upper bounds

intcon = 1;   %  x(1) must be an integer, number of firms

warning('off', 'MATLAB:nearlySingularMatrix')

%%%%%% settings for genetic algorithm

A   = []; b   = [];
Aeq = []; beq = []; nonlcon = [];

addpath(genpath('/scratch/jab772/compecom')); 
matlabpool open 8
%matlabpool open 2

%gaoptions = gaoptimset(gaoptions,'UseParallel','always','Display','iter','InitialPopulation',xstart);
gaoptions = gaoptimset('UseParallel', 'always', 'TimeLimit',86400,'Display','iter','InitialPopulation',x);
x = ga(@objective, nvars, A, b, Aeq, beq, lb, ub, [], intcon, gaoptions);

%options = gaoptimset('Generations',100,'Vectorized','on','UseParallel','always');    
%[x,fval,exitflag,output] = ga(@(x)objective(x,params),nvars,A,b,Aeq,beq,lb,ub,nonlcon,intcon);

xga=x;

save Results xga

matlabpool close

x     = x';
xga = x;   % save in case fminsearch goes wild

x = fminsearch('objective', x, optimset('Display', 'iter')); 

save Results xga x

