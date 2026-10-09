function [S_m,K_m,D_m] = Sim_CN_DNN_1(DNN,par_solid,dim_model)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
H = par_solid(1); %cm Altura de la pila
T = par_solid(2); %días, tiempo de simulación
theta_ini = par_solid(3); %Humedad inicial (cm3/cm3)
theta_s = par_solid(4); %Humedad de saturación (cm3/cm3)
theta_r = par_solid(5); %Humedad resifual (cm3/cm3)
Ks = par_solid(6); % (cm/día) Conductividad hidráulica de saturación

Dt = dim_model(1); %1/día Dt es por cada hora
Dz = dim_model(2); %cm

%Cuerpo
D_H = 24; %conversión día/hora
A = 308; %área de irrigación m2

Nz = H/Dz; %número de puntos de simulación en eje z
Nt = T/Dt; %número de puntos de simulación en eje t
Nz = round(Nz);
par_solid = [theta_s;theta_r;Ks]; %DATA LAB
dim_model = [Nz; Dz; Dt];


%Vector de caudales diarios (un caudal por cada día)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ESTA OPCION SE USA PARA CAUDALES IGUALES DURANTE LOS 44 DIAS
%q0_init = 25; %m3/d Caudal de ingreso inicial (VAR DE DECISION)
%Q_0_vect = ones(T,1)*q0_init; %44 caudales 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ESTA OPCION SE USA PARA UN CAUDAL POR DIA
Q_0_vect = [12;19;19;0;0;0;12;32;36;39;46;27;23;24;27;30;29;28;28;27;36;35;26;26;26;25;22;26;24;24;26;22;22;25;27;27;23;19;20;19;16;17;21;30]; %44 caudales
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Vector de caudales de salida
qout = zeros(length(Q_0_vect),1);

%Redes Neuronales
net_k_d = DNN{1};
max_K = DNN{2};
max_D = DNN{3};
par_key=[max_D,max_K];



%Matriz de resultados
theta_m_dnn = zeros(Nz+1,Nt+1);
S_m = zeros(Nz+1,Nt+1);
K_m = zeros(Nz+1,Nt+1);
D_m = zeros(Nz+1,Nt+1);
%Vector de humedades iniciales
theta_v = ones(Nz+1,1)*theta_ini;
theta_act = theta_v;
k = 1;
theta_m_dnn(:,k) = theta_act;
S_act = S_theta(theta_act,theta_r,theta_s,Nz);
[Y_1_2] = predict(net_k_d,S_act);
K_ini = Y_1_2(:,1).*max_K;
D_ini = Y_1_2(:,2).*max_D;

S_m(:,k) = S_act;
K_m(:,k) = K_ini;
D_m(:,k) = D_ini;


n_d = length(Q_0_vect);

T_dias = zeros(length(T),1);


for d = 1:n_d %paso del tiempo en días
    
    T_dias(d,1) = d;
    
    q0 = Q_0_vect(d,1); %m3/d Caudal de ingreso (VAR DE DECISION) se cambia el caudal en cada día (cada vez que pasa el ciclo for superior
    
    R = (q0/A)*100; % tasa de riego en (cm/d)
    
    for h = 1:D_H %contador de horas por cada día
        
        k = (d-1)*D_H + h; %contador de horas correlativas
        
        %A(theta_act)*theta_p = b(theta_act) + c(theta_act) (primera iter) -->
        %predictor
        
        %A(theta_p)*theta_c = b(theta_act) + c(theta_p) (pimera iter) -->
        
        b = Pred_correct_b_dnn(theta_act,R,par_solid,dim_model,net_k_d,par_key);
        theta_v = theta_act;
        N = 20;
        theta_v_p_iter = zeros(Nz+1,N);
        
        for j = 1:N
            theta_v_p_iter(:,j) = Pred_correct_dnn(theta_v,b,R,par_solid,dim_model,net_k_d,par_key);
            theta_v = theta_v_p_iter(:,j);
        end
        S_v = S_theta(theta_v,theta_r,theta_s,Nz);
        [Y_1_2] = predict(net_k_d,S_v);
        K = Y_1_2(:,1).*max_K;
        D = Y_1_2(:,1).*max_D;

        S_m(:,k+1) = S_v;
        K_m(:,k+1) = K;
        D_m(:,k+1) = D;

        theta_act = theta_v;
        theta_m_dnn(:,k+1) = theta_act;
        
        %%%%%CALCULO DE CAUDAL DE SALIDA EN CADA h%%%%
        %S = S_theta(theta_act,theta_r,theta_s,Nz);
        %K = K_theta(S,Ks,n,Nz,method);
        %qout(1,k) = K(Nz+1,1);
        
        
    end

    qout(d,1) = K(Nz+1,1)*(A/100);% m3/dia
end

end