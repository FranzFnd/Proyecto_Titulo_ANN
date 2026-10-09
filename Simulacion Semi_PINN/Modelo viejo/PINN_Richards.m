clc
clearvars
%Parámetros
H = 120; %cm Altura de la pila
T = 44; %días, tiempo de simulación
theta_ini = 0.14; %Humedad inicial (cm3/cm3)
theta_s = 0.33; %Humedad de saturación (cm3/cm3)
theta_r = 0; %Humedad resifual (cm3/cm3)
Ks = 170; % (cm/día) Conductividad hidráulica de saturación
alpha = 0.035; %(1/cm) Parámetro VGM
n = 2.267; %Parámetro VGM
method = 2; %método VGM
Dt = 1/24; %1/día Dt es por cada hora
Dz = 15; %cm

D_H = 24; %conversión día/hora
A = 308; %área de irrigación m2

Q_0_vect = [12;19;19;0;0;0;12;32;36;39;46;27;23;24;27;30;29;28;28;27;36;35;26;26;26;25;22;26;24;24;26;22;22;25;27;27;23;19;20;19;16;17;21;30];

Z_vect = [0:Dz:H];
T_vect = [0:Dt:T];
Nz = H/Dz; %número de puntos de simulación en eje z
Nt = T/Dt; %número de puntos de simulación en eje t

par_solid = [theta_s;theta_r;Ks]; %DATA LAB
par_model = [alpha;n;method]; 
dim_model = [Nz; Dz; Dt];


%Vector de caudales de salida
qout = zeros(length(Q_0_vect),1);

theta_v = ones(Nz+1,1)*theta_ini;
theta_act = theta_v;
k = 1;
theta_m_pred(:,k) = theta_act;
input = theta_m_pred;
%Extraccion de Datos de entrenamiento
load("Datos theta.mat");
output = theta_m;

%Definir arquitectura de NN
%%%%%%%NN K y D
numBlocks = 1;
fcOutputSize = 11;

fcBlock = [
    fullyConnectedLayer(fcOutputSize)
    tanhLayer];

layers = [
    featureInputLayer(1)
    repmat(fcBlock,[numBlocks 1])
    fullyConnectedLayer(2)];

net_k_d = dlnetwork(layers);

net_k_d = dlupdate(@double,net_k_d);


function [loss,gradients] = modelLoss(net_k_d,input,theta_m)
    %Parámetros
    H = 120; %cm Altura de la pila
    T = 44; %días, tiempo de simulación
    theta_ini = 0.14; %Humedad inicial (cm3/cm3)
    theta_s = 0.33; %Humedad de saturación (cm3/cm3)
    theta_r = 0; %Humedad resifual (cm3/cm3)
    Ks = 170; % (cm/día) Conductividad hidráulica de saturación
    alpha = 0.035; %(1/cm) Parámetro VGM
    n = 2.267; %Parámetro VGM
    method = 2; %método VGM
    Dt = 1/24; %1/día Dt es por cada hora
    Dz = 15; %cm

    D_H = 24; %conversión día/hora
    A = 308; %área de irrigación m2

    Q_0_vect = [12;19;19;0;0;0;12;32;36;39;46;27;23;24;27;30;29;28;28;27;36;35;26;26;26;25;22;26;24;24;26;22;22;25;27;27;23;19;20;19;16;17;21;30];

    Z_vect = [0:Dz:H];
    T_vect = [0:Dt:T];
    Nz = H/Dz; %número de puntos de simulación en eje z
    Nt = T/Dt; %número de puntos de simulación en eje t

    par_solid = [theta_s;theta_r;Ks]; %DATA LAB
    par_model = [alpha;n;method]; 
    dim_model = [Nz; Dz; Dt];
    % Paso 1: Inferir y1 y y2 con las redes neuronales
    %procesar datos antes y despues

    [y_1_2] = forward(net_k_d,input);
    

    theta_out = ModelRichard_PINN(K,D,par_solid,dim_model,Q_0_vect);

    loss1 = l2loss(theta_out,theta_m);
    %loss2 = l2loss(qout,qout0);

    loss = loss1;


    % Calculate gradients with respect to the learnable parameters.
    gradients = dlgradient(loss,net_k_d.Learnables);

end

%Se definene la funcion de optimizacion, iteraciones maximas, tolerancia y
%paso
solverState = lbfgsState;
maxIterations = 1500;
gradientTolerance = 1e-7;
stepTolerance = 1e-7;
%Se ordenan los datos en array
theta_input = dlarray(input',"CB");

theta_obj = dlarray(output,"TCB");

accfun = dlaccelerate(@modelLoss);
lossFcn = @(net_k_d) dlfeval(accfun,net_k_d,theta_input,theta_obj);
%Display
monitor = trainingProgressMonitor( ...
    Metrics="TrainingLoss", ...
    Info=["Iteration" "GradientsNorm" "StepNorm"], ...
    XLabel="Iteration");
%Entrenamiento de la red neuronal
iteration = 0;
while iteration < maxIterations && ~monitor.Stop
    iteration = iteration + 1;
    [net_k_d, solverState] = lbfgsupdate(net_k_d,lossFcn,solverState);

    updateInfo(monitor, ...
        Iteration=iteration, ...
        GradientsNorm=solverState.GradientsNorm, ...
        StepNorm=solverState.StepNorm);

    recordMetrics(monitor,iteration,TrainingLoss=solverState.Loss);

    monitor.Progress = 100*iteration/maxIterations;

    if solverState.GradientsNorm < gradientTolerance || ...
            solverState.StepNorm < stepTolerance || ...
            solverState.LineSearchStatus == "failed"
        break
    end

end



