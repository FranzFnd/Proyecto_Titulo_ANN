function [theta_m1, D_m1, K_m1] = theta_m1_theta_PINN(theta_0,D0,K0,Dz,R,theta_r,theta_s,net)
%THETA_M1_THETA Summary of this function goes here
%   Detailed explanation goes here


theta_m1 = theta_0 + (R - K0)*Dz/(D0); %despeje de theta_-1 de CB superior
if theta_m1 >= theta_s
    diff = theta_m1-theta_s;
    theta_m1 = theta_m1 - 1.05*diff;
end
D_theta_dif = theta_s - theta_r;
S = (theta_m1 - theta_r)/D_theta_dif;

[K_m1,D_m1] = forward(net,theta_m1);


end

