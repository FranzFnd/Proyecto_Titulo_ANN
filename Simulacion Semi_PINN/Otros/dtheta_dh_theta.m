function [dtheta_dh] = dtheta_dh_theta(h,theta_s,theta_r,n,alpha,Nz,method)
%DTHETA_DH_THETA Summary of this function goes here
%   Detailed explanation goes here

dtheta_dh = zeros(1,1);
D_theta = theta_s - theta_r;

if method == 1 %método VG
    m = 1 - 1/n;
else %método VGM
    m = 1 - 2/n;
end
mult = D_theta*n*alpha*(-m);
%mult = -m*n*(alpha^n)*D_theta;

for i = 1
    %dtheta_dh(i) = mult*((abs(h(i))^n)*(1/h(i))*(1/((1+(alpha^n)*abs(h(i))^n)^(m+1))));
    dtheta_dh(i) = mult*((alpha*h(i))^(n-1))/((1+(alpha*h(i))^n)^(m+1));
end





end

