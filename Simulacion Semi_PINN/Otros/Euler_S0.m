function [Se0] = Euler_S0(Se, D, K, Rj, Param)


dt=Param(1);dtau=Param(2);H=Param(3);dchi=Param(4);dtheta=Param(5);

Se = extractdata(Se);
Se = Se';
D = (D).*1000;
K = (K).*10;

i=1;

Se0 = (((Rj-K(i))*H*dchi)/(D(i)*dtheta))+Se(i+1);


end