function [S_m,K_m,D_m] = Sim_CN_VG(par_solid,par_vg,dim_model)
%Sim_CN_VG utiliza el modelo de van genuchten para obtener los valores de
%K, h, y D a partir de un theta dado y crank-nicholson para simular la
%ecuacion diferencial parcial y obtener el theta del momento siguiente
%Parametros

H = par_solid(1); %cm Altura de la pila
T = par_solid(2); %días, tiempo de simulación
theta_ini = par_solid(3); %Humedad inicial (cm3/cm3)
theta_s = par_solid(4); %Humedad de saturación (cm3/cm3)
theta_r = par_solid(5); %Humedad resifual (cm3/cm3)
Ks = par_solid(6); % (cm/día) Conductividad hidráulica de saturación
alpha = par_vg(1); %(1/cm) Parámetro VGM
n = par_vg(2); %Parámetro VGM
method = par_vg(3); %método VGM
Dt = dim_model(1); %1/día Dt es por cada hora
Dz = dim_model(2); %cm

%Cuerpo 

D_H = 24; %conversión día/hora
A = 308; %área de irrigación m2


Nz = H/Dz; %número de puntos de simulación en eje z
Nt = T/Dt; %número de puntos de simulación en eje t

par_solid = [theta_s;theta_r;Ks]; %DATA LAB
par_model = [alpha;n;method]; 
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


%Matriz de resultados
theta_m = zeros(Nz+1,Nt+1);

S_m = zeros(Nz+1,Nt+1);
K_m = zeros(Nz+1,Nt+1);
dtheta_dh_theta_m = zeros(Nz+1,Nt+1);
D_m = zeros(Nz+1,Nt+1);

%Vector de humedades iniciales
theta_v = ones(Nz+1,1)*theta_ini;
theta_act = theta_v;
k = 1;
theta_m(:,k) = theta_act;

S = S_theta(theta_v,theta_r,theta_s,Nz);
ache = h_theta(S,n,alpha,Nz,method);
K = K_theta(S,Ks,n,Nz,method);
dtheta_dh = dtheta_dh_theta(ache,theta_s,theta_r,n,alpha,Nz,method);
D = D_theta(K,dtheta_dh,Nz);

S_m(:,k) = S;
K_m(:,k) = K;
dtheta_dh_theta_m(:,k) = dtheta_dh;
D_m(:,k) = D;

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
        
        b = Pred_correct_b(theta_act,R,par_solid, par_model,dim_model);
        theta_v = theta_act;
        N = 20;
        theta_v_p_iter = zeros(Nz+1,N);
        
        for j = 1:N
            theta_v_p_iter(:,j) = Pred_correct(theta_v,b,R,par_solid, par_model,dim_model);
            theta_v = theta_v_p_iter(:,j);
        end
        S = S_theta(theta_v,theta_r,theta_s,Nz);
        ache = h_theta(S,n,alpha,Nz,method);
        K = K_theta(S,Ks,n,Nz,method);
        dtheta_dh = dtheta_dh_theta(ache,theta_s,theta_r,n,alpha,Nz,method);
        D = D_theta(K,dtheta_dh,Nz);

        S_m(:,k+1) = S;
        K_m(:,k+1) = K;
        dtheta_dh_theta_m(:,k+1) = dtheta_dh;
        D_m(:,k+1) = D;
        
        theta_act = theta_v;
        theta_m(:,k+1) = theta_act;
        
        %%%%%CALCULO DE CAUDAL DE SALIDA EN CADA h%%%%
        %S = S_theta(theta_act,theta_r,theta_s,Nz);
        %K = K_theta(S,Ks,n,Nz,method);
        %qout(1,k) = K(Nz+1,1);
        
        
    end
    
    S = S_theta(theta_act,theta_r,theta_s,Nz);
    K = K_theta(S,Ks,n,Nz,method);% cm/dia
    qout(d,1) = K(Nz+1,1)*(A/100);% m3/dia
end
end