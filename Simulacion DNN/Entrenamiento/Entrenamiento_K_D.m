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

%Generacion de datos
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

Se_val = Se_val';
D_val = D_val';
K_val = K_val';

%isequal(sort([trainIdx,ValIdx]), 1:numel(train_G_Idx))



% Inicializar redes neuronales, se genera la arquitectura de la red, una
% vez creada se puede verificar con analyzeNetwork(net_k_d)
numBlocks = 2;
fcOutputSize = 10;

fcBlock = [
            fullyConnectedLayer(fcOutputSize)
            tanhLayer];

        layers = [
            featureInputLayer(1)
            repmat(fcBlock,[numBlocks 1])
            fullyConnectedLayer(2)];
% functionLayer(@(x) log(1+exp(x)),Description="softplusLayer") para usar
% softplus al final de la arquitectura
        net_k_d = dlnetwork(layers);

        net_k_d = dlupdate(@double,net_k_d);

function [loss,gradients] = modelLoss(net_k_d,S,K0,D0)

    % Make predictions with the initial conditions.
    S_ini = S;
    [Y_1_2] = forward(net_k_d,S_ini);
    
    K_out = Y_1_2(1,:);
    D_out = Y_1_2(2,:);
    
    loss1 = l2loss(K_out,K0);
    loss2 = l2loss(D_out,D0);
    
    
    loss = loss1+loss2;
    
    
    % Calculate gradients with respect to the learnable parameters.
    gradients = dlgradient(loss,net_k_d.Learnables);

end

%Se definene la funcion de optimizacion, iteraciones maximas, toleeancia y
%paso
solverState = lbfgsState;
maxIterations = 1500;
gradientTolerance = 1e-5;
stepTolerance = 1e-6;
%Se ordenan los datos en array
Se = dlarray(Se_train,"CB");
D0 = dlarray(D_train,"CB");
K0 = dlarray(K_train,"CB");
trainLossHistory = zeros(1, maxIterations);
valLossHistory = zeros(1, maxIterations);

accfun = dlaccelerate(@modelLoss);
lossFcn = @(net_k_d) dlfeval(accfun,net_k_d,Se,D0,K0);
%Display
monitor = trainingProgressMonitor( ...
    Metrics="TrainingLoss", ...
    Info=["Iteration" "GradientsNorm" "StepNorm" "Loss"], ...
    XLabel="Iteration");
%Entrenamiento de la red neuronal
iteration = 0;
while iteration < maxIterations && ~monitor.Stop
    iteration = iteration + 1;
    [net_k_d, solverState] = lbfgsupdate(net_k_d,lossFcn,solverState);

    trainLossHistory(iteration) = solverState.Loss;

        % Forward pass and loss computation for validation
    Y_1_2 = predict(net_k_d, Se_val);
    K_val_pred = Y_1_2(:,1);
    D_val_pred = Y_1_2(:,2);
    valLoss = mean((K_val - K_val_pred).^2) + mean((D_val - D_val_pred).^2); % Mean Squared Error
    valLossHistory(iteration) = valLoss;

    updateInfo(monitor, ...
        Iteration=iteration, ...
        GradientsNorm=solverState.GradientsNorm, ...
        StepNorm=solverState.StepNorm, ...
        Loss = solverState.Loss);

    recordMetrics(monitor,iteration,TrainingLoss=solverState.Loss);

    monitor.Progress = 100*iteration/maxIterations;

        if solverState.GradientsNorm < gradientTolerance || ...
            solverState.StepNorm < stepTolerance || ...
            solverState.LineSearchStatus == "failed"
        break
        end


end
disp("Optimización completada.");
DNN_K_D = cell(1,3);
DNN_K_D{1} = net_k_d;
DNN_K_D{2} = max_K;
DNN_K_D{3} = max_D;
save("DNN_K_D","DNN_K_D")

%Validacion

figure(1);
plot(1:iteration, trainLossHistory(1:iteration), '-b', 'LineWidth', 2); hold on;
plot(1:iteration, valLossHistory(1:iteration), '-r', 'LineWidth', 2);
xlabel('época');
ylabel('Pérdida');
title('Desempeño de entrenamiento y validación');
legend('Pérdida Entrenamiento', 'Pérdida Validación');
grid on;

Se_train = Se_train';
Y_1_2 = predict(net_k_d, Se_train);
K_test_pred = Y_1_2(:,1);
D_test_pred = Y_1_2(:,2);
K_test_pred = K_test_pred';
D_test_pred = D_test_pred';

%K
SSres = sum((K_train - K_test_pred).^2);
SStot = sum((K_train - mean(K_train)).^2);
R2_K = 1 - (SSres / SStot);

figure(2);
scatter(K_train, K_test_pred, 'filled');
hold on;
plot(K_train, K_train, 'r', 'LineWidth', 2); 
xlabel('Valores reales');
ylabel('Valores predichos');
title(sprintf('Predicciones de entrenamiento vs Real K (R²: %.4f)', R2_K));
grid on;

errors_K = K_train - K_test_pred;
figure(3);
histogram(errors_K, 'Normalization', 'probability');
xlabel('Error de predicción');
ylabel('Probablidad');
title('Histograma de errores de K');
grid on;


%D
SSres = sum((D_train - D_test_pred).^2);
SStot = sum((D_train - mean(D_train)).^2);
R2_D = 1 - (SSres / SStot);

figure(4);
scatter(D_train, D_test_pred, 'filled');
hold on;
plot(D_train, D_train, 'r', 'LineWidth', 2); 
xlabel('Valores reales');
ylabel('Valores predichos');
title(sprintf('Predicciones de entrenamiento vs Real D (R²: %.4f)', R2_D));
grid on;

errors_D = D_train - D_test_pred;
figure(5);
histogram(errors_D, 'Normalization', 'probability');
xlabel('Error de predicción');
ylabel('Probablidad');
title('Histograma de errores de D');
grid on;

%-------------------------------------------
%Validacion datos

Se_val = Se_val;
Y_1_2 = predict(net_k_d, Se_val);
K_test_pred = Y_1_2(:,1);
D_test_pred = Y_1_2(:,2);
K_test_pred = K_test_pred;
D_test_pred = D_test_pred;


%K
SSres = sum((K_val - K_test_pred).^2);
SStot = sum((K_val - mean(K_val)).^2);
R2_K_val = 1 - (SSres / SStot);

    figure(6);
scatter(K_val, K_test_pred, 'filled');
hold on;
plot(K_val, K_val, 'r', 'LineWidth', 2); % Ideal fit line
xlabel('Valores reales');
ylabel('Valores predichos');
title(sprintf('Predicciones de validación vs Real K (R²: %.4f)', R2_K_val));
grid on;

errors_K = K_val - K_test_pred;
    figure(7);
histogram(errors_K, 'Normalization', 'probability');
xlabel('Error de predicción');
ylabel('Probablidad');
title('Histograma de errores de K');
grid on;


%D
SSres = sum((D_val - D_test_pred).^2);
SStot = sum((D_val - mean(D_val)).^2);
R2_D_val = 1 - (SSres / SStot);

    figure(8);
scatter(D_val, D_test_pred, 'filled');
hold on;
plot(D_val, D_val, 'r', 'LineWidth', 2); % Ideal fit line
xlabel('Valores reales');
ylabel('Valores predichos');
title(sprintf('Predicciones de validación vs Real D (R²: %.4f)', R2_D_val));
grid on;

errors_D = D_val - D_test_pred;
    figure(9);
histogram(errors_D, 'Normalization', 'probability');
xlabel('Error de predicción');
ylabel('Probablidad');
title('Histograma de errores de D');
grid on;