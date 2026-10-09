

clc; clear; close all;

% Parámetros iniciales para usar la funcion toy y para generar los datos de
% prueba de salida
%N_z = 100; % Número total de puntos espaciales
%N_p = 20; % Número de puntos con mediciones
%z = linspace(0, 1, N_z); % Discretización espacial
%indices_p = randperm(N_z, N_p); % Índices de posiciones medidas
%t = 0:0.1:5; % Tiempo de simulación
%u = 1; % Entrada constante como ejemplo
%Parámetros
H = 550; %cm Altura de la pila
T = 44; %días, tiempo de simulación
theta_ini = 0.14; %Humedad inicial (cm3/cm3)
theta_s = 0.33; %Humedad de saturación (cm3/cm3)
theta_r = 0; %Humedad resifual (cm3/cm3)
dtheta = theta_s-theta_r;
Ks = 170; % (cm/día) Conductividad hidráulica de saturación
Dt = 1/24; %1/día Dt es por cada hora
Dz = 10; %cm
Dchi = 10;
Dtau = 1/24;

D_H = 24; %conversión día/hora
A = 308; %área de irrigación m2

Z_vect = [0:Dz:H];
T_vect = [0:Dt:T];
Nz = H/Dz; %número de puntos de simulación en eje z
Nt = T/Dt; %número de puntos de simulación en eje t
Nz = round(Nz);
Param = [Dt, Dtau, H, Dchi, dtheta];
% q0_init = [25];%cm/dia
% Q_0_vect = ones(T,1)*q0_init;
Q_0_vect = [12;19;19;0;0;0;12;32;36;39;46;27;23;24;27;30;29;28;28;27;36;35;26;26;26;25;22;26;24;24;26;22;22;25;27;27;23;19;20;19;16;17;21;30];

%Simualcion de S_m K_m D_m
[S_m,K_m,D_m] = Sim_CN_VG(par_solid,par_vg,dim_model);
% Inicializar redes neuronales, se genera la arquitectura de la red, una
% vez creada se puede verificar con analyzeNetwork(net_k_d)
numBlocks = 2;
fcOutputSize = 20;

fcBlock = [
            fullyConnectedLayer(fcOutputSize)
            tanhLayer];

        layers = [
            featureInputLayer(1)
            repmat(fcBlock,[numBlocks 1])
            fullyConnectedLayer(2)
            softmaxLayer];
% functionLayer(@(x) log(1+exp(x)),Description="softplusLayer") para usar
% softplus al final de la arquitectura
        net_k_d = dlnetwork(layers);

        net_k_d = dlupdate(@double,net_k_d);



%Cálculo de Pérdida y Gradientes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [loss, gradients] = modelLoss(net_k_d, Se_m_full, Nz, Nt, Q_0_vect, Param)

    %primer instante de x_simulado
    Se_sim = zeros(Nz+1,Nt+1);%se genera una matriz
    % Se_sim = dlarray(Se_sim,'CB');%se transforma a dlarray object
    Se_sim(:, 1) = Se_m_full(:, 1); % Condición inicial

    n_d = length(Q_0_vect);
    D_H = 24; %conversión día/hora
    A = 308; %área de irrigación m2
    
    for d = 1:n_d %paso del tiempo en días
    
        q0 = Q_0_vect(d,1); %m3/d Caudal de ingreso (VAR DE DECISION) se cambia el caudal en cada día (cada vez que pasa el ciclo for superior
    
        Rj = (q0/A)*100; % tasa de riego en (cm/d)
    
        for h = 1:D_H %contador de horas por cada día
        
            k = (d-1)*D_H + h; %contador de horas correlativas   
            if k == 1
                   
            else
                Se_prev = Se_sim(:, k-1);%se calcula el x del momento k con el x del momeno k-1 osea el x previo
                Se_prev = dlarray(Se_prev','CB');
                Y_1_2 = forward(net_k_d,Se_prev);
                D = Y_1_2(1,:);
                D = extractdata(D);
                D = D';
                K = Y_1_2(2,:);
                K = extractdata(K);
                K = K';
                Se_sim(:, k) = euler_Forward(Se_prev, D, K, Rj, Param);%aqui va la funcion del modelo
            end
        end
    end
    Se_sim=dlarray(Se_sim,'CB');
    loss = l2loss(Se_sim,Se_m_full); % Distancia cuadrática
    loss = extractdata(loss);
    % Extraer parámetros de la red como tabla
    Learnables = net_k_d.Learnables;
    
    % Obtener los parámetros de la tabla
    Layer = Learnables{:, 1};
    Parameter = Learnables{:, 2};
    
    numElementos = numel(Learnables{:, 3});
    Valordlarray = Learnables{:, 3}; % Datos en formato dlarray dentro de celdas
    
    % Convertir valores dlarray a tipo double
    Valores = cellfun(@extractdata, Valordlarray, 'UniformOutput', false);
    
    % Aplicar gradiente y guardar cada resultado por separado
    Values = cell(numElementos, 1);
    for i = 1:numElementos
        Matriz = Valores{i};  % Extraer matriz numérica
        gradiente = gradient(Matriz, loss);
        Values{i, 1} = gradiente;
    
    end

    % Crear tabla final con gradientes aplicados
    gradients = table(Layer, Parameter, Values, ...
                           'VariableNames', {'Layer', 'Parameter', 'Value'});
    %gradients = dlgradient(loss,net_k_d.Learnables);
    

end



%Se definene la funcion de optimizacion, iteraciones maximas, tolerancia y
%paso
solverState = lbfgsState;
maxIterations = 1500;
gradientTolerance = 1e-7;
stepTolerance = 1e-7;
%Se ordenan los datos en array

accfun = dlaccelerate(@modelLoss);
lossFcn = @(net_k_d) dlfeval(@modelLoss,net_k_d, Se_m_full, Nz, Nt, Q_0_vect, Param);
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

    updateInfo(monitor, ...
        Iteration=iteration, ...
        GradientsNorm=solverState.GradientsNorm, ...
        StepNorm=solverState.StepNorm,Loss=solverState.Loss);

    recordMetrics(monitor,iteration,TrainingLoss=solverState.Loss);

    monitor.Progress = 100*iteration/maxIterations;

    if solverState.GradientsNorm < gradientTolerance || ...
            solverState.StepNorm < stepTolerance || ...
            solverState.LineSearchStatus == "failed"
        disp("Optimización Fallida.");
        break
    end

end


disp("Optimización completada.");
save("DNN_K_D","net_k_d")




