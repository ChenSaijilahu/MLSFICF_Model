clc;clear all;
%alpha = 1.0;   % ��
lambda = 0.2;   % �� = 0.0, 0.1, 0.2, 0.3
gamma = 0.1;    % �� = 0.0, 0.1, 0.2, 0.3

P0 = 0.90;      %�� = 1, 0.9, 0.95, 0.85
P1 = 0.95;      %�� = 1, 0.9, 0.95, 0.85
P2 = 0.95;      %�� = 1, 0.9, 0.95, 0.85

V1 = 6.75;
V2 = 7.91;
C1 = 0.13;
C2 = 1.57;
lc = 5.0;
v_max = 2;

syms h n
Optimal_Velocity_F = @(h) V1+V2*tanh(C1*(h-lc)-C2);
Optimal_Velocity_B = @(h) -(V1+V2*tanh(C1*(h-lc)-C2));
Optimal_Velocity_F_Diff = matlabFunction(diff(Optimal_Velocity_F(h)));   
Optimal_Velocity_B_Diff = matlabFunction(diff(Optimal_Velocity_B(h))); 
Optimal_Velocity_F_Diff3 = matlabFunction(diff(Optimal_Velocity_F(h),h,3));   
Optimal_Velocity_B_Diff3 = matlabFunction(diff(Optimal_Velocity_B(h),h,3)); 

cc1 = P1*(6/7*1+6/49*3+1/49*5)-(1-P1);
cc2 = P2*(6/7*1+6/49*2+1/49*3)-(1-P2);
cc3 = P1*(6/7*1+6/49*7+1/49*19)+(1-P1);
cc4 = P2*(6/7*1+6/49*4+1/49*9)+(1-P2);

dd1 = @(h) P0*Optimal_Velocity_F_Diff(h)+(1-P0)*Optimal_Velocity_B_Diff(h);
dd2 = @(h) P0*(6/7*1+6/49*3+1/49*5)*Optimal_Velocity_F_Diff(h)-(1-P0)*Optimal_Velocity_B_Diff(h);
dd3 = @(h) P0*(6/7*1+6/49*7+1/49*19)*Optimal_Velocity_F_Diff(h)+(1-P0)*Optimal_Velocity_B_Diff(h);
dd4 = @(h) P0*Optimal_Velocity_F_Diff3(h)+(1-P0)*Optimal_Velocity_B_Diff3(h);
dd5 = @(h) P0*(6/7*1+6/49*15+1/49*65)*Optimal_Velocity_F_Diff(h)-(1-P0)*Optimal_Velocity_B_Diff(h);
dd6 = @(h) P0*(6/7*1+6/49*3+1/49*5)*Optimal_Velocity_F_Diff3(h)-(1-P0)*Optimal_Velocity_B_Diff3(h);

H0 = 17.07692308;
alpha_x = @(h) 2*((1-gamma)*dd1(h)*dd1(h)-lambda*dd1(h))/dd2(h);

alpha_c = alpha_x(H0)

g1 = @(h) 1/6*dd3(h) + 1/2*lambda*dd1(h)/alpha_c*cc1 + gamma*dd1(h)*dd1(h)/alpha_c*cc2;
g2 = @(h) -1/6*dd4(h);
g3 = @(h) (1-gamma)*dd1(h)*dd1(h)/alpha_c - lambda*dd1(h)/alpha_c;
g4 = @(h) (2*dd1(h)*(1-gamma)-lambda)/alpha_c*(1/6*dd3(h)+1/2*lambda*dd1(h)*cc1/alpha_c+gamma*dd1(h)*dd1(h)*cc2/alpha_c)-(1/24*dd5(h)+1/6*lambda*dd1(h)*cc3/alpha_c+1/2*gamma*dd1(h)*dd1(h)*cc4/alpha_c);
g5 = @(h) 1/12*(2*(2*dd1(h)*(1-gamma)-lambda)/alpha_c*dd4(h)-dd6(h));

c = 5*g2(H0)*g3(H0) / (2*g2(H0)*g4(H0) - 3*g1(H0)*g5(H0))
% G = [g1(H0),g2(H0),g3(H0),g4(H0),g5(H0),c,alpha_c]


for H = 0.1:0.1:36.0
    Num = round(H*10);
    Y_MLSFICF3(Num) = alpha_x(H);
    YY_MLSFICF3(Num) = alpha_c/(1+(H-H0)^2*g2(H0)/(g1(H0)*c));
end

figure(1);
X = 0.1:0.1:36.0;
plot(X,Y_MLSFICF3,'b+');hold on;
plot(X,YY_MLSFICF3,'r+');hold off;

axis([0 35 0 2.5]);