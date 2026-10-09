function [Se_next] = RK4_Richards(Se, D, K, Rj, Param)
%UNTITLED3 Summary of this function goes here
%   PARAMETRS
    dt=Param(1);dtau=Param(2);H=Param(3);dchi=Param(4);dtheta=Param(5);
    Se_next = zeros(length(Se),1);
    N = length(Se);
    K1 = zeros(1,N);%Inicializa los vectores para guardar los valores de K1, K2, K3 Y K4
    K2 = zeros(1,N);
    K3 = zeros(1,N);
    K4 = zeros(1,N);
    
    Se0 = Se(2) + (H*dchi)*(-K(1)+Rj)/(dtheta*D(1)); %Calculo de condicion inicial cuando i=1

    i=1;
    K1(1) = (dt/(H^2))*(((D(i+1)-D(i))/(2*dchi))*((Se(i+1)-Se0)/(2*dchi)) + D(i)*((Se(i+1)-2*Se(i)+Se0)/(dchi^2)) -((K(i+1)-K(i))/(2*dchi))*(H/dtheta));
    D_temp = D + dtau/2*K1(1);
    K_temp = K + dtau/2*K1(1);
    Se0_temp = Se0 + dtau/2*K1(1); 
    Se_temp = Se + dtau/2*K1(1);
    K2(1) = (dt/(H^2))*(((D_temp(i+1)-D_temp(i))/(2*dchi))*((Se_temp(i+1)-Se0_temp)/(2*dchi)) + D_temp(i)*((Se_temp(i+1)-2*Se_temp(i)+Se0_temp)/(dchi^2)) -((K_temp(i+1)-K_temp(i))/(2*dchi))*(H/dtheta));
    
    D_temp = D + dtau/2*K2(1);
    K_temp = K + dtau/2*K2(1);
    Se0_temp = Se0 + dtau/2*K2(1); 
    Se_temp = Se + dtau/2*K2(1);
    K3(1) = (dt/(H^2))*(((D_temp(i+1)-D_temp(i))/(2*dchi))*((Se_temp(i+1)-Se0_temp)/(2*dchi)) + D_temp(i)*((Se_temp(i+1)-2*Se_temp(i)+Se0_temp)/(dchi^2)) -((K_temp(i+1)-K_temp(i))/(2*dchi))*(H/dtheta));
   
    D_temp = D + dtau*K3(1);
    K_temp = K + dtau*K3(1);
    Se0_temp = Se0 + dtau*K3(1); 
    Se_temp = Se + dtau*K3(1);
    K4(1) = (dt/(H^2))*(((D_temp(i+1)-D_temp(i))/(2*dchi))*((Se_temp(i+1)-Se0_temp)/(2*dchi)) + D_temp(i)*((Se_temp(i+1)-2*Se_temp(i)+Se0_temp)/(dchi^2)) -((K_temp(i+1)-K_temp(i))/(2*dchi))*(H/dtheta));

    Se_next(1) = Se(1) + (dtau/6)*(K1(1)+2*K2(1)+2*K3(1)+K4(1));

    %RK4 para la simulacion
    for i = 2:N-1
    K1(i) = (dt/(H^2))*(((D(i+1)-D(i-1))/(2*dchi))*((Se(i+1)-Se(i-1))/(2*dchi)) + D(i)*((Se(i+1)-2*Se(i)+Se(i-1))/(dchi^2)) -((K(i+1)-K(i-1))/(2*dchi))*(H/dtheta));
    end
    D_temp = D + dtau/2*K1;
    K_temp = K + dtau/2*K1;
    Se_temp = Se + dtau/2*K1;
    for i = 2:N-1
    K2(i) = (dt/(H^2))*(((D_temp(i+1)-D_temp(i-1))/(2*dchi))*((Se_temp(i+1)-Se_temp(i-1))/(2*dchi)) + D_temp(i)*((Se_temp(i+1)-2*Se_temp(i)+Se_temp(i-1))/(dchi^2)) -((K_temp(i+1)-K_temp(i))/(2*dchi))*(H/dtheta));
    end
    D_temp = D + dtau/2*K2;
    K_temp = K + dtau/2*K2;
    Se_temp = Se + dtau/2*K2;
    for i = 2:N-1
    K3(i) = (dt/(H^2))*(((D_temp(i+1)-D_temp(i-1))/(2*dchi))*((Se_temp(i+1)-Se_temp(i-1))/(2*dchi)) + D_temp(i)*((Se_temp(i+1)-2*Se_temp(i)+Se_temp(i-1))/(dchi^2)) -((K_temp(i+1)-K_temp(i))/(2*dchi))*(H/dtheta));
    end
    D_temp = D + dtau*K3;
    K_temp = K + dtau*K3;
    Se_temp = Se + dtau*K3;
    for i = 2:N-1
    K4(i) = (dt/(H^2))*(((D_temp(i+1)-D_temp(i-1))/(2*dchi))*((Se_temp(i+1)-Se_temp(i-1))/(2*dchi)) + D_temp(i)*((Se_temp(i+1)-2*Se_temp(i)+Se_temp(i-1))/(dchi^2)) -((K_temp(i+1)-K_temp(i))/(2*dchi))*(H/dtheta));
    end
    for i = 2:N-1
    Se_next(i,1) = Se(i) + (dtau/6)*(K1(i)+2*K2(i)+2*K3(i)+K4(i));
    end
    Se_next(N) = Se_next(N-1);%Condicion final
end