function [S] = S_theta(theta,theta_r,theta_s,Nz)
%S_THETA 
% S = (theta-theta_r)/(theta_s - theta_r)
D_theta = theta_s - theta_r;

S = (1/D_theta)*(theta - theta_r*ones(Nz+1,1));
end

