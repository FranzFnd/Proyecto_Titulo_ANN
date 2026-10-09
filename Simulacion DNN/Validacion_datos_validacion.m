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

[Nz,Nt] = size(S_m);
Se_in = S_m(1,:);
D_in = D_m(1,:);
K_in = K_m(1,:);

max_D =  max(D_in);
D_in = D_in./max_D;

max_K =  max(K_in);
K_in = K_in./max_K;

train_Idx = randperm(numel(Se_in), round(numel(Se_in)*.7));
Val_Idx = find(~ismember(1:numel(Se_in), train_Idx));

Se_train = Se_in(train_Idx);
D_train = D_in(train_Idx);
K_train = K_in(train_Idx);

Se_val = Se_in(Val_Idx);
D_val = D_in(Val_Idx);
K_val = K_in(Val_Idx);

load('DNN_K');
load('DNN_D.mat');
net_k = DNN_K{1};
max_K2 = DNN_K{2};
net_d = DNN_D{1};
max_D2 = DNN_D{2};

%Datos train

%Validacion de los modelos
%K


Se_train = Se_train';
K_test_pred = predict(net_k, Se_train);
D_test_pred = predict(net_d, Se_train);
K_test_pred = K_test_pred';
D_test_pred = D_test_pred';

%K
SSres = sum((K_train - K_test_pred).^2);
SStot = sum((K_train - mean(K_train)).^2);
R2_K = 1 - (SSres / SStot);

    figure(3);
scatter(K_train, K_test_pred, 'filled');
hold on;
plot(K_train, K_train, 'r', 'LineWidth', 2); % linea de ajuste ideal
xlabel('Valores reales');
ylabel('Valores predichos');
title(sprintf('Predicciones de entrenamiento vs Real K (R²: %.4f)', R2_K));
grid on;

errors_K = K_train - K_test_pred;
    figure(4);
histogram(errors_K, 'Normalization', 'probability');
xlabel('Error de predicción');
ylabel('Probablidad');
title('Histograma de errores de K');
grid on;


%D
SSres = sum((D_train - D_test_pred).^2);
SStot = sum((D_train - mean(D_train)).^2);
R2_D = 1 - (SSres / SStot);

    figure(5);
scatter(D_train, D_test_pred, 'filled');
hold on;
plot(D_train, D_train, 'r', 'LineWidth', 2); 
xlabel('Valores reales');
ylabel('Valores predichos');
title(sprintf('Predicciones de entrenamiento vs Real D (R²: %.4f)', R2_D));
grid on;

errors_D = D_train - D_test_pred;
    figure(6);
histogram(errors_D, 'Normalization', 'probability');
xlabel('Error de predicción');
ylabel('Probablidad');
title('Histograma de errores de D');
grid on;

%________________________________________________________________________
%Datos validacion
Se_val = Se_val';
K_test_pred = predict(net_k, Se_val);
D_test_pred = predict(net_d, Se_val);
K_test_pred = K_test_pred';
D_test_pred = D_test_pred';

%K
SSres = sum((K_val - K_test_pred).^2);
SStot = sum((K_val - mean(K_val)).^2);
R2_K = 1 - (SSres / SStot);

    figure(7);
scatter(K_val, K_test_pred, 'filled');
hold on;
plot(K_val, K_val, 'r', 'LineWidth', 2); % Ideal fit line
xlabel('Valores reales');
ylabel('Valores predichos');
title(sprintf('Predicciones de validación vs Real K (R²: %.4f)', R2_K));
grid on;

errors_K = K_val - K_test_pred;
    figure(8);
histogram(errors_K, 'Normalization', 'probability');
xlabel('Error de predicción');
ylabel('Probablidad');
title('Histograma de errores de K');
grid on;


%D
SSres = sum((D_val - D_test_pred).^2);
SStot = sum((D_val - mean(D_val)).^2);
R2_D = 1 - (SSres / SStot);

    figure(9);
scatter(D_val, D_test_pred, 'filled');
hold on;
plot(D_val, D_val, 'r', 'LineWidth', 2); % Ideal fit line
xlabel('Valores reales');
ylabel('Valores predichos');
title(sprintf('Predicciones de validación vs Real D (R²: %.4f)', R2_D));
grid on;

errors_D = D_val - D_test_pred;
    figure(10);
histogram(errors_D, 'Normalization', 'probability');
xlabel('Error de predicción');
ylabel('Probablidad');
title('Histograma de errores de D');
grid on;