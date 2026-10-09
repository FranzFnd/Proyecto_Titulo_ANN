function [theta_m1, D_m1, K_m1] = theta_m1_theta(theta_0,D0,K0,Dz,R,n,alpha,Ks,theta_r,theta_s,method)
%THETA_M1_THETA Summary of this function goes here
%   Detailed explanation goes here

theta_m1 = theta_0 + (R - K0)*Dz/D0; %despeje de theta_-1 de CB superior
D_theta_dif = theta_s - theta_r;
S = (theta_m1 - theta_r)/D_theta_dif;

K_0 = K_theta(S,Ks,n,0,method);
h_0 = h_theta(S,n,alpha,0,method);
dtheta_dh_0 = dtheta_dh_theta(h_0,theta_s,theta_r,n,alpha,0,method);

D_m1 = D_theta(K_0,dtheta_dh_0,0);
K_m1 = K_0;

end

