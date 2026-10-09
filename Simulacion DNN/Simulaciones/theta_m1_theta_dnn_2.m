function [theta_m1, D_m1, K_m1] = theta_m1_theta_dnn_2(theta_0,D0,K0,Dz,R,theta_r,theta_s,net_k,net_D,par_key)
%THETA_M1_THETA Summary of this function goes here
%   Detailed explanation goes here
max_D = par_key(1); max_K = par_key(2);

theta_m1 = theta_0 + (R - K0)*Dz/(D0); %despeje de theta_-1 de CB superior
if theta_m1 >= theta_s
    diff = theta_m1-theta_s;
    theta_m1 = theta_m1 - 1.05*diff;
end
D_theta_dif = theta_s - theta_r;
S = (theta_m1 - theta_r)/D_theta_dif;

K_0 = predict(net_k,S).*max_K;

D_m1 = predict(net_D,S).*max_D;
K_m1 = K_0;


end