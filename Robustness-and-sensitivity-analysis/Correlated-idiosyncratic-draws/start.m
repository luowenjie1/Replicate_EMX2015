clear;
clc;
rng('default');

% Assigned Parameters

b1 = -0.503;       % slope coefficient and intercept in regression of inverse markups on market shares
b0 = 0.637;

S  = 1000;         % 1/2 number of sectors
L  = 1;            % inelastic labor supply in Home
Ls = 1;            % inelastic labor supply in Foreign

N  = 279;          % maximum number of producers per sector


% Calibrated Parameters

gamma = 10.4971041171121;                   % within-sector elasticity of substitution   
theta = (1/gamma - b1/b0*(1-1/gamma))^(-1); % across-sector elasticity of substitution
xi_x  = 4.58079215954791;                   % Pareto shape, idiosyncratic productivity
xi_z  = 0.505927117530769;                  % Pareto shape, sector productivity      
zeta  = 0.0428446336934024;                 % geometric parameter, number producers per sector
FD    = 0.0044022548579538;                 % fixed cost of domestic operations
FX    = 0.06;                               % fixed cost of export operations
tau   = 0.1285;                             % net trade cost 
rho_z = 0.934;                              % Kendall correlation, *sectoral productivity*

rho_x = 0.05;                               % Kendall correlation, *idiosyncratic productivity*


tausave  = tau;   % save since need later
   
if 1

fprintf('\n');    
display('Computing Autarky Equilibrium')    
fprintf('\n');

tau = 100000;
comp_trade_elasticity = 0;  % numerically compute trade elasticity if =1

equilibrium;
        
Asave = A; 
Alosssave = log(Aeff/A)*100;
    
save autarkysaved Asave Alosssave

% Uncomment the following to report more statistics

%model_moments;
%markup_moments;
     
%break 
        
end

fprintf('\n');
display('Computing Equilibrium with Trade')
fprintf('\n');

tau = tausave; 
comp_trade_elasticity = 1;
        
equilibrium
        
load autarkysaved  % CAREFUL: IF PARAMETERS CHANGE, THIS NEEDS TO BE UPDATED 

fprintf('A increase (gains from trade), *100      = %7.3f \n', log(A/Asave)*100);  
fprintf('Due to markups                           = %7.3f \n', Alosssave - log(Aeff/A)*100);  
fprintf('Trade elasticity                         = %7.3f \n', sigma);  
fprintf('Import share                             = %7.3f \n', impshare);
fprintf('Fraction exporters                       = %7.3f \n', agg_fexporters);  

% Uncomment the following to report more statistics

%model_moments;
%markup_moments;
%domestic_and_import_markups;         
     

        

