clear all
clc

%%%%% parameters to calibrate

N       = 110;                  % number of firms that draw
ebar    = 1.4866;          % productivity of type 2
su      = 0.1955 ;         % standard deviation of Home/Foreign productivity gap
tau     = 1.1331; 

om1     = [ 0.1216, 0.0283 , 0.3388]; 
om2     = [ 0.3271, 0.1262, 0.2030]; 

pom     = [0.4775, 0.8366,  0.1584,... 
           0.2276, 0.2994, 0.7672,... 
            0.9008, 0.7324, 0.0385];
       
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

