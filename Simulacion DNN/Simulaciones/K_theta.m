function [K] = K_theta(S,Ks,n,Nz,method)
%K_THETA Summary of this function goes here
%   Detailed explanation goes here

K = zeros(Nz+1,1);

if method == 1 %método VG
    m = 1 - 1/n;
    for i = 1 : Nz+1
        K(i) = Ks*sqrt(S(i))*(1-(1-S(i)^(1/m))^m)^2;
    end
else %método VGM
    m = 1 - 2/n;
    delta = 3 + 2/(m*n);
    for i = 1 : Nz+1
        K(i) = Ks*(S(i))^delta;
    end
end


end

