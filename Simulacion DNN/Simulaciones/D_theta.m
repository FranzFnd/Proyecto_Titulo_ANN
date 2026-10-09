function [D] = D_theta(K,dtheta_dh,Nz)
%D_THETA Summary of this function goes here
%   Detailed explanation goes here

D = zeros(Nz+1,1);

for i = 1:Nz+1
    D(i) = -K(i)/dtheta_dh(i);
end




end

