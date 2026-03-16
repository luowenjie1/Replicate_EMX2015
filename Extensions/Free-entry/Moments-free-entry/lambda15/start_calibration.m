clear all
clc

%%%%% parameters to calibrate

N       = 140;                                  % number of firms that draw
ebar    = 1.7477;                          % productivity of type 2
su      = 0.1706;                         % standard deviation of Home/Foreign productivity gap
tau     = 1.1269; 

om1     = [0.1085 , 0.0128,  0.3006]; 
om2     = [0.0642, 0.0976, 0.0516]; 

pom     = [ 0.9147, 0.5660,  0.2033, ...
           0.4803,0.4149 ,0.4887 , ...
           0.5512, 0.6701, 0.4749];

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

