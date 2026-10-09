clc; clear; close all;
%Creacion de red neuronal para K

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
percentage_below_0_55 = (sum(Se_in < 0.55) / length(Se_in)) * 100;
% Define the edges of the intervals
edges = 0:0.05:1;

% Calculate the histogram
[counts, ~] = histcounts(Se_in, edges);

% Calculate the percentage of total data points in each interval
percentage_counts = (counts / length(Se_in)) * 100;

% Create the bar graph
figure;
bar(edges(1:end-1), percentage_counts, 'histc');
xlabel('Intervalos de Se');
ylabel('Porcentaje de datos en el intervalo (%)');
title('Histograma de distribución de valores de Se');
xlim([0 1]);
grid on;