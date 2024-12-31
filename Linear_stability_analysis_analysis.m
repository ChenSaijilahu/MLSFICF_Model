clc; clear all;

alpha = 1.0;      % α
lambda = 0.2;     % λ = 0.0, 0.1, 0.2, 0.3
gamma = 0.1;      % γ = 0.0, 0.1, 0.2, 0.3

P = 0.9;  % η = 1,0.9,0.95,0.85
P1 = 0.9; % μ
P2 = 0.9; % ρ
q = 7;

V1 = 6.75;
V2 = 7.91;
C1 = 0.13;
C2 = 1.57;
lc = 5.0;

syms h n
Optimal_Velocity_OV = @(h) V1+V2*tanh(C1*(h-lc)-C2);
Optimal_Velocity_FVD = @(h) V1+V2*tanh(C1*(h-lc)-C2);
Optimal_Velocity_FVDA = @(h) V1+V2*tanh(C1*(h-lc)-C2);
Optimal_Velocity_BLVD_F = @(h) V1+V2*tanh(C1*(h-lc)-C2);
Optimal_Velocity_BLVD_B = @(h) -(V1+V2*tanh(C1*(h-lc)-C2));
Optimal_Velocity_MLSFICF1_F = @(h) V1+V2*tanh(C1*(h-lc)-C2);
Optimal_Velocity_MLSFICF1_B = @(h) -(V1+V2*tanh(C1*(h-lc)-C2));
Optimal_Velocity_MLSFICF3_F = @(h) V1+V2*tanh(C1*(h-lc)-C2);
Optimal_Velocity_MLSFICF3_B = @(h) -(V1+V2*tanh(C1*(h-lc)-C2));

Optimal_Velocity_OV_Diff = matlabFunction(diff(Optimal_Velocity_OV(h)));   
Optimal_Velocity_FVD_Diff = matlabFunction(diff(Optimal_Velocity_FVD(h))); 
Optimal_Velocity_FVDA_Diff = matlabFunction(diff(Optimal_Velocity_FVDA(h))); 
Optimal_Velocity_BLVD_Diff_F = matlabFunction(diff(Optimal_Velocity_BLVD_F(h))); 
Optimal_Velocity_BLVD_Diff_B = matlabFunction(diff(Optimal_Velocity_BLVD_B(h))); 
Optimal_Velocity_MLSFICF_1_Diff_F = matlabFunction(diff(Optimal_Velocity_MLSFICF1_F(h))); 
Optimal_Velocity_MLSFICF_1_Diff_B = matlabFunction(diff(Optimal_Velocity_MLSFICF1_B(h))); 
Optimal_Velocity_MLSFICF_3_Diff_F = matlabFunction(diff(Optimal_Velocity_MLSFICF3_F(h))); 
Optimal_Velocity_MLSFICF_3_Diff_B = matlabFunction(diff(Optimal_Velocity_MLSFICF3_B(h))); 

D1_BLVD = @(h) P*Optimal_Velocity_BLVD_Diff_F(h) + (1-P)*Optimal_Velocity_BLVD_Diff_B(h);
D2_BLVD = @(h) P*Optimal_Velocity_BLVD_Diff_F(h) - (1-P)*Optimal_Velocity_BLVD_Diff_B(h);
D1_1 = @(h) P*Optimal_Velocity_MLSFICF_1_Diff_F(h) + (1-P)*Optimal_Velocity_MLSFICF_1_Diff_B(h);
D2_1 = @(h) P*Optimal_Velocity_MLSFICF_1_Diff_F(h) - (1-P)*Optimal_Velocity_MLSFICF_1_Diff_B(h);
D1_3 = @(h) P*Optimal_Velocity_MLSFICF_3_Diff_F(h) + (1-P)*Optimal_Velocity_MLSFICF_3_Diff_B(h);
D2_3 = @(h) P*(6/7*1+6/49*3+1/49*5)*Optimal_Velocity_MLSFICF_3_Diff_F(h) - (1-P)*Optimal_Velocity_MLSFICF_3_Diff_B(h);

for H = 0.1:0.1:36.0
    Num = round(H*10);
    Y_OV(Num) = 2*Optimal_Velocity_OV_Diff(H);
    Y_FVD(Num) = 2*(Optimal_Velocity_FVD_Diff(H)-lambda);
    Y_FVDA(Num) = 2*(Optimal_Velocity_FVDA_Diff(H)*(1-gamma)-lambda);
    Y_BLVD(Num) = 2*D1_BLVD(H)*(D1_BLVD(H) - lambda)/D2_BLVD(H);
    Y_MLSFICF_1(Num) = 2*(D1_1(H)*D1_1(H) - lambda*D1_1(H) - gamma*D1_1(H)*D1_1(H))/D2_1(H);
    Y_MLSFICF_3(Num) = 2*(D1_3(H)*D1_3(H) - lambda*D1_3(H) - gamma*D1_3(H)*D1_3(H))/D2_3(H);
end

X = 0.1:0.1:36.0;
plot(X,Y_OV,'y+');hold on;
plot(X,Y_FVD,'m*');hold on;
plot(X,Y_FVDA,'c-');hold on;
plot(X,Y_BLVD,'r-');hold on;
plot(X,Y_MLSFICF_1,'g:');hold on;
plot(X,Y_MLSFICF_3,'b');hold off;
axis([0 36 0 2.5]);
