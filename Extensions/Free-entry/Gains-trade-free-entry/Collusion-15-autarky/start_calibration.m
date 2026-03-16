clear all
clc

%%%%% parameters to calibrate

N       = 150;                  % number of firms that draw
om1     = [0.02, 0.0988, 0.35]; 
om2     = [0.02, 0.0620, 0.12]; 
ebar    = 2.0;                  % productivity of type 2
su      = 0.25;                 % standard deviation of Home/Foreign productivity gap
tau     = 1.15; 

pom     = 1/9*ones(1,9); 

x      = [N, ebar, su, tau, om1, om2,pom];

nvars  = length(x);

%%%%% constraints on x

lb     = [100,   1.1,  0.05,  1.01,...
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
matlabpool open 12
%matlabpool open 2

%gaoptions = gaoptimset(gaoptions,'UseParallel','always','Display','iter','InitialPopulation',xstart);
gaoptions = gaoptimset('UseParallel', 'always','StallGenLimit', 3600*2, 'Display','iter','InitialPopulation',x);
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

