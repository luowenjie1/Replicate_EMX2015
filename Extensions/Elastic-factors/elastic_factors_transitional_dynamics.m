function F = elastic_labor_dynamics2(X,parameters,Abar,Mbar,Kbar,K0,Cbar,Lbar)

alpha = parameters(1) ; %% output elasticity of capital in aggregate production function
delta = parameters(2) ; %% capital depreciation
beta  = parameters(3) ; %% discount factor
chi   = parameters(4) ; %% weight on leisure (chosen so L=1)
eta   = parameters(5) ; %% inverse Frisch elasticity labor supply


Cpath = [X(:,1);Cbar];
Kpath = [X(:,2);Kbar];

T = length(X(:,1));

Lpath = (((1-alpha)/chi)*(Abar/Mbar)*(Kpath.^(alpha))./Cpath).^(1/(eta+alpha)); %% closed-form solution for L given C,K

%%%%% calculate equation errors

F = zeros(T,2);

for t=1:T,
   
    F(t,1) = Abar*(Kpath(t)^alpha)*(Lpath(t)^(1-alpha)) + (1-delta)*Kpath(t) - Cpath(t) -Kpath(t+1); %% resource constraint
    
    F(t,2) = (Cpath(t+1)/Cpath(t)) - beta*(alpha*(Abar/Mbar)*(Kpath(t+1)^(alpha-1))*(Lpath(t+1)^(1-alpha))+1-delta); %% euler equation    

end


%norm(F);
     
    




