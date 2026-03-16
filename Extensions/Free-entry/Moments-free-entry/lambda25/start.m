clear all;
clc;

global gamma theta tau ebar N S

% Assigned Parameters

gamma   = 10.4971;                    
theta   = 1.2351;

S       = 500;                  % number of industries: must be even
L       = 1;                    % inelastic labor supply in Home
Ls      = 1;                    % inelastic labor supply in Foreign

lambda  = 0.25;                 % prob. all large guys collude
pom     = 1;                    % ok for now


% Calibrated

N       = 160;                             % number of firms that draw
ebar    = 1.7176;                          % productivity of type 2
su      = 0.1737;                          % standard deviation of Home/Foreign productivity gap
tau     = 1.1356; 

om1     = [0.1583 ; 0.0116; 0.3131]; 
om2     = [0.1535; 0.0587;0.2596]; 

pom     = [0.4681; 0.3000; 0.4239; ...
           0.5701;0.6364 ;0.6787 ; ...
           0.4236; 0.7156; 0.6289];
       
pom     = pom/sum(pom); 

comp_armington = 1;                 % compute armington elasticity

[uu, wu] = qnwunif(7, 0, 1);
uu       = exp(norminv(uu, 0, su)); 

equilibrium
model_moments

break
tau=1000;
comp_armington = 0;                 % compute armington elasticity
equilibrium

% break
% % equilibrium_check
% % break
% % equilibrium
% % break
% % % Choose what exactly you want to report
% 
% switch 'equilibrium'
% 
%     case 'equilibrium'
%          
%         if 0
%         
%         tau = 100000; 
%         tausave = tau; 
%         
%         equilibrium;
%         
%         Asave = A; 
%         Alosssave = log(Aeff/A)*100;
%     
%         save autarkysaved Asave Alosssave
%             
%         end
%         
%         %tau = 0.4045;       % 20% import share
%         %tausave = tau; 
%         
%         %tau = 0.129;        % Taiwan's import share
%         %tausave = tau; 
%         
%         %tau = 0.7058; 
%         %tausave = tau; 
%         
%         tau = 0; 
%         tausave = tau; 
%         
%         equilibrium
%         
%         load autarkysaved
%         
%         fprintf('A increase (gains from trade), *100      = %7.3f \n', log(A/Asave)*100);  
%         fprintf('Due to markups                           = %7.3f \n', Alosssave - log(Aeff/A)*100);  
%         fprintf('Trade elasticity                         = %7.3f \n', - Arm);  
%         fprintf('Import share                             = %7.3f \n', impshare);  
% 
%         
%     case 'model_moments'
%         
%         equilibrium;
%         
%         fprintf('\n');
%         fprintf('\n');
%                 
%         model_moments;
% 
%         
% y = nanmean(abs(100*loss));
% disp('value of objective')
% disp(y);       
%  
%     case 'data_moments' % only run this if you want to re-calculate data moments
%         
%         data_moments;
%      
%         
% end
%         
%     
