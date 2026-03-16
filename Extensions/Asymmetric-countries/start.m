clear all;
clc;
rng('default');

% Assigned Parameters

b1    = -0.503;       % slope coefficient and intercept in regression of inverse markups on market shares
b0    = 0.637;

S     = 1000;         % 1/2 number of sectors
L     = 1;            % inelastic labor supply in Home
Ls    = 1;            % inelastic labor supply in Foreign

Abars = 10;           % relative productivity  in Foreign 

N     = 279;          % maximum number of producers per sector


% Calibrated Parameters
 
gamma = 10.25;                              % within-sector elasticity of substitution   
theta = (1/gamma - b1/b0*(1-1/gamma))^(-1); % across-sector elasticity of substitution
xi_x  = 4.58079215954791;                   % Pareto shape, idiosyncratic productivity
xi_z  = 0.505927117530769;                  % Pareto shape, sector productivity      
zeta  = 0.0428446336934024;                 % geometric parameter, number producers per sector
FD    = 0.0044022548579538;                 % fixed cost of domestic operations
FX    = 0.30;                               % fixed cost of export operations
tau   = 1.66;                               % net trade cost 
rho   = 0.59975;                            % Kendall correlation for Gumbel copula 

FDs   = FD;

tausave  = tau;   % save since need later
 

% Asymmetric country experiments have benchmark parameters except that for

%Ls = 2         : tau = 0.245, rho = 0.94425             
%Ls = 10        : tau = 0.50 , rho = 0.955225             
%Abars = 2      : tau = 0.32 , rho = 0.86    , fx = 0.25, gamma = 10.25
%Abars = 10     : tau = 1.66 , rho = 0.59975 , fx = 0.30, gamma = 10.25

% Benchmark values, for reference:

%gamma = 10.4971041171121;                    
%FX    = 0.2025;                             
%tau   = 0.1285;                            
%rho   = 0.935325;                            



if 1

fprintf('\n');    
display('Computing Autarky Equilibrium')    
fprintf('\n');

tau = 100000;
comp_trade_elasticity = 0;  % numerically compute trade elasticity if =1
        
equilibrium_autarky;
        
Asave = A;
Assave = As;
        
Aesave = Ae;
Asesave = Ase;
        
Alosssave = log(Ae/A)*100;
Aslosssave = log(Ase/As)*100;
        
save autarkysaved Asave Assave Aesave Asesave Alosssave Aslosssave
        
%break 
end

fprintf('\n');
display('Computing Equilibrium with Trade')
fprintf('\n');

tau = tausave; 
comp_trade_elasticity = 1;
        
equilibrium
        
load autarkysaved  % CAREFUL: IF PARAMETERS CHANGE, THIS NEEDS TO BE UPDATED 

fprintf('A increase (gains from trade), *100 (H)      = %7.3f \n', log(A/Asave)*100); 
fprintf('A increase (gains from trade), *100 (F)      = %7.3f \n', log(As/Assave)*100); 
fprintf('\n');
fprintf('Due to markups                      (H)      = %7.3f \n', Alosssave - log(Ae/A)*100); 
fprintf('Due to markups                      (F)      = %7.3f \n', Aslosssave - log(Ase/As)*100);
fprintf('\n');
fprintf('Trade elasticity                    (H)      = %7.3f \n', sigma);  
fprintf('Trade elasticity                    (F)      = %7.3f \n', sigmas);
fprintf('\n');
fprintf('Import share                        (H)      = %7.3f \n', impshare); 
fprintf('Import share                        (F)      = %7.3f \n', impshares); 
fprintf('\n');
fprintf('Fraction exporters                  (H)      = %7.3f \n', agg_fexporters);
fprintf('Fraction exporters                  (F)      = %7.3f \n', agg_fexporterss); 

% Uncomment the following to report more statistics

%model_moments_asymmetric;
%markup_moments_asymmetric;         
     


