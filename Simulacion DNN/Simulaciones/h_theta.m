function [h] = h_theta(S,n,alpha,Nz,method)
%H_THETA Summary of this function goes here
%   Detailed explanation goes here

if method == 1 %método VG
    m = 1 - 1/n;
else %método VGM
    m = 1 - 2/n;
end

h = zeros(Nz+1,1);
alpha_inv = 1/alpha;

for i = 1:Nz+1
    h(i) =  alpha_inv*((1/S(i))^(1/m) - 1)^(1/n);
end

end

