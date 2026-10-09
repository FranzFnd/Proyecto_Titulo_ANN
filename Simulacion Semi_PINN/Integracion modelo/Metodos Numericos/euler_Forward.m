function [Se_next] = euler_Forward(Se, D, K, Rj, Param)
    
    dt=Param(1);dtau=Param(2);H=Param(3);dchi=Param(4);dtheta=Param(5);
    % alpha = 0.035; %(1/cm) Parámetro VGM
    % n = 2.267; %Parámetro VGM
    % method = 2; %método VGM
    % Ks = 170; % (cm/día) Conductividad hidráulica de saturación
    % theta_s = 0.33; %Humedad de saturación (cm3/cm3)
    % theta_r = 0; %Humedad resifual (cm3/cm3)
    % dtheta = theta_s-theta_r;
    N = length(Se);
    %Codigo necesario para usar en funcion Training_net_D_K.m
    % Se0 = extractdata(Se0);
    % Se = extractdata(Se);
    % Se = Se';
    % D = (D).*1000;
    % K = (K).*10;

   
    % Inicializa el nuevo vector Se(:, j+1)
    Se_next = zeros(length(Se),1);
    for i = 1:N
        if i == 1 %Condicion de borde i=1
            
            Se0 = Se(2) + (H*dchi)*(-K(1)+Rj)/(dtheta*D(1));


            % h0 = h_theta(Se0,n,alpha,0,method);
            % K0 = K_theta(Se0,Ks,n,0,method);
            % dtheta_dh = dtheta_dh_theta(h0,theta_s,theta_r,n,alpha,0,method);
            % D0 = D_theta(K,dtheta_dh,0);    

            Term1 = ((D(i+1)-D(i))/(dchi))*((Se(i+1)-Se0)/(2*dchi));
            Term2 = D(i)*((Se(i+1)-2*Se(i)+Se0)/(dchi^2));
            Term3 = ((K(i+1)-K(i))/(dchi))*(H/dtheta);
            Se_next(i,1)= Se(i)+(dt*dtau/(H^2))*(Term1+Term2-Term3);
        elseif i == N %condicion de borde i=N

            Se_next(i,1) = Se(i);
        else %Desde i=2:N-1
            Term1 = ((D(i+1)-D(i-1))/(2*dchi))*((Se(i+1)-Se(i-1))/(2*dchi));
            Term2 = D(i)*((Se(i+1)-2*Se(i)+Se(i-1))/(dchi^2));
            Term3 = ((K(i+1)-K(i-1))/(2*dchi))*(H/dtheta);
            Se_next(i,1)= Se(i)+(dt*dtau/(H^2))*(Term1+Term2-Term3);
        end
    end



end