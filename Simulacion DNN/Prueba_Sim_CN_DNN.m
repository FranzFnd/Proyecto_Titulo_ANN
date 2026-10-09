clc; clear; close all;

H = 550; %cm Altura de la pila
T = 44; %días, tiempo de simulación
theta_ini = 0.14; %Humedad inicial (cm3/cm3)
theta_s = 0.33; %Humedad de saturación (cm3/cm3)
theta_r = 0; %Humedad resifual (cm3/cm3)
dtheta = theta_s-theta_r;
Ks = 170; % (cm/día) Conductividad hidráulica de saturación
alpha = 0.035; %(1/cm) Parámetro VGM
n = 2.267; %Parámetro VGM
method = 2; %método VGM
Dt = 1/24; %1/día Dt es por cada hora
Dz = 10; %cm

par_solid = [H,T,theta_ini,theta_s,theta_r,Ks];
par_vg = [alpha,n,method];
dim_model = [Dt,Dz];

Z_vect = 0:Dz:H;
T_vect = 0:Dt:T;

[S_m,K_m,D_m] = Sim_CN_VG(par_solid,par_vg,dim_model);
load('DNN_K_D.mat');
load('DNN_K');
load('DNN_D.mat');
net_k_d = DNN_K_D{1};
max_K1 = DNN_K_D{2};
max_D1 = DNN_K_D{3};
net_k = DNN_K{1};
max_K2 = DNN_K{2};
net_d = DNN_D{1};
max_D2 = DNN_D{2};
[S_m1,K_m1,D_m1] = Sim_CN_DNN_1(DNN_K_D,par_solid,dim_model);
[S_m2,K_m2,D_m2] = Sim_CN_DNN_2(DNN_K,DNN_D,par_solid,dim_model);

    figure(1)
hold on
mesh(T_vect, Z_vect, S_m)
xlabel('tiempo (d)')
ylabel('profundidad (m)')
zlabel('humedad (-)')
title('Evolución del perfil de humedad para R dado')


figure(2)
hold on
mesh(T_vect, Z_vect, S_m1)
xlabel('tiempo (d)')
ylabel('profundidad (m)')
zlabel('humedad (-)')
title('Evolución del perfil de humedad para R dado')



figure(3)
hold on
mesh(T_vect, Z_vect, S_m2)
xlabel('tiempo (d)')
ylabel('profundidad (m)')
zlabel('humedad (-)')
title('Evolución del perfil de humedad para R dado')

%Para resultados de red con K y D juntos
SSres = sum((S_m - S_m1).^2, "all");
SStot = sum((S_m - mean(S_m1)).^2,"all");
R2_S_1 = 1 - (SSres / SStot);

S_m_vector = S_m(:);
S_m1_vector = S_m1(:);

    figure(4);
scatter(S_m_vector, S_m1_vector, 'filled');
hold on;
plot(S_m_vector, S_m_vector, 'r', 'LineWidth', 2); % Ideal fit line
xlabel('Valores actuales');
ylabel('Valores predichos');
title(sprintf('Predicciones Se vs VGM Se (R²: %.4f)', R2_S_1));
grid on;

%Para resultados de red con K y D separados
SSres = sum((S_m - S_m2).^2, "all");
SStot = sum((S_m - mean(S_m2)).^2,"all");
R2_S_2 = 1 - (SSres / SStot);

S_m_vector = S_m(:);
S_m2_vector = S_m2(:);

    figure(5);
scatter(S_m_vector, S_m2_vector, 'filled');
hold on;
plot(S_m_vector, S_m_vector, 'r', 'LineWidth', 2); % Ideal fit line
xlabel('Valores actuales');
ylabel('Valores predichos');
title(sprintf('Predicciones Se vs VGM Se (R²: %.4f)', R2_S_2));
grid on;