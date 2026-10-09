% File: pdae_nn_ad.m

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
Dchi = Dz/H;
Dtau = (Dt-0)/(T-0);

D_H = 24; %conversión día/hora
A = 308; %área de irrigación m2

Z_vect = [0:Dz:H];
T_vect = [0:Dt:T];
Nz = H/Dz; %número de puntos de simulación en eje z
Nt = T/Dt; %número de puntos de simulación en eje t
Nz = round(Nz);
Param = [Dt, Dtau, H, Dchi, dtheta]';
% q0_init = [25];%cm/dia
% Q_0_vect = ones(T,1)*q0_init;
Q_0_vect = [12;19;19;0;0;0;12;32;36;39;46;27;23;24;27;30;29;28;28;27;36;35;26;26;26;25;22;26;24;24;26;22;22;25;27;27;23;19;20;19;16;17;21;30];

load('Datos S.mat');
Se_m_full=S_m;
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

                Se_prev = Se_sim(:, k);%se calcula el x del momento k con el x del momeno k-1 osea el x previo
                Se_prev = dlarray(Se_prev','CB');
                Y_1_2 = forward(net_k_d,Se_prev);
                D = Y_1_2(1,:);
                D = extractdata(D);
                D = D';
                K = Y_1_2(2,:);
                K = extractdata(K);
                K = K';
                % Usar codigo 103 hasta 111 para obtener Se0 junto con D0 Y
                % K0 de la red neuronal
                % Se0 = Euler_S0(Se_prev, D, K, Rj, Param);
                % Se0 = dlarray(Se0,'CB');
                % Y_1_2 = forward(net_k_d,Se0);
                % D0 = Y_1_2(1,:);
                % D0 = extractdata(D0);
                % D0 = D0';
                % K0 = Y_1_2(2,:);
                % K0 = extractdata(K0);
                % K0 = K0';
                Se_sim(:, k+1) = euler_Forward(Se_prev, D, K, Rj, Param);%aqui va la funcion del modelo
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
     % Finite Difference Gradient Approximation

end



% Optimización, este tema es aparte 
numEpochs = 1000;
numIterationsPerEpoch = 10;
iteration = 0;
epoch = 0;
averageGrad = [];
averageSqGrad = [];
numIterations = numEpochs * numIterationsPerEpoch;
monitor = trainingProgressMonitor(Metrics="Loss",Info="Epoch",XLabel="Epoch");
while epoch < numEpochs && ~monitor.Stop
    epoch = epoch + 1;
    i = 0;
    disp(["epoca numero", epoch]);
    while i < numIterationsPerEpoch && ~monitor.Stop
        i = i + 1;
        iteration = iteration + 1;

        % Evaluate the model loss and gradients using dlfeval and the
        % modelLoss function.
        [loss,gradients] = dlfeval(@modelLoss,net_k_d, Se_m_full, Nz, Nt, Q_0_vect, Param);

        % Update the network parameters using the Adam optimizer.
        [net_k_d,averageGrad,averageSqGrad] = adamupdate(net_k_d,gradients,averageGrad,averageSqGrad,iteration);

        % Update the training progress monitor.
        recordMetrics(monitor,iteration,Loss=loss);
        updateInfo(monitor,Epoch=epoch + " of " + numEpochs);
        monitor.Progress = 100 * iteration/numIterations;
    end
end

%Validacion


% Visualización

disp("Optimización completada.");
save("DNN_K_D","net_k_d")




