function [Theta_next] = Euler_no_adim(theta, D, K, Rn, Param)

dt=Param(1);dz=Param(2);
N = length(theta);
Theta_next = zeros(length(theta),1);
    for i = 1:N-1
        if i==1
        theta0 = ((Rn-K(i))*2*dz)/D(i)+theta(i+1);
        Term1 = ((D(i+1)-D(1))/(2*dz))*((theta(i+1)-theta0)/(2*dz));
        Term2 = D(i)*((theta(i+1)-2*theta(i)+theta0)/(dz^2));
        Term3 = ((K(i+1)-K(1))/(2*dz));
        Theta_next(i,1)= theta(i)+dt*(Term1+Term2-Term3);
        else
        Term1 = ((D(i+1)-D(i-1))/(2*dz))*((theta(i+1)-theta(i-1))/(2*dz));
        Term2 = D(i)*((theta(i+1)-2*theta(i)+theta(i-1))/(dz^2));
        Term3 = ((K(i+1)-K(i-1))/(2*dz));
        Theta_next(i,1)= theta(i)+dt*(Term1+Term2-Term3);
        end
    end
    i=N;
    Theta_next(i,1)= theta(i);
    
end