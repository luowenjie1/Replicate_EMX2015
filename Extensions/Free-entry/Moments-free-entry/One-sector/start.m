clear all;
clc;

global gamma theta tau ebar N S

% Assigned Parameters

gamma   = 10.4971;                    
theta   = 1.2351;
S       = 1000;                 % number of industries: must be even
L       = 1;                    % inelastic labor supply in Home
Ls      = 1;                    % inelastic labor supply in Foreign

lambda  = 0.0;                  % prob. all large guys collude
pom     = 1;                    % ok for now


% Calibrated

N       = 176;                  % number of firms that draw
om1     = 0.0972513900; 
om2     = 0.1734328984; 
ebar    = 1.838675631;          % productivity of type 2
su      = 0.1469919780;         % standard deviation of Home/Foreign productivity gap
tau     = 1.1339970969; 


comp_armington = 1;             % compute armington elasticity

[uu, wu] = qnwunif(7, 0, 1);
uu       = exp(norminv(uu, 0, su)); 


equilibrium
model_moments

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
