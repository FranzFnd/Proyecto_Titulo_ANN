function [theta_out] = ModelRichard_PINN(K_matrx,D_matrx,par_solid,dim_model,Q_0_vect)
    
    theta_s = par_solid(1); theta_r = par_solid(2); Ks = par_solid(3);
    Nz = dim_model(1); Dz = dim_model(2); Dt = dim_model(3);
    %Parámetros
    H = 120; %cm Altura de la pila
    T = 44; %días, tiempo de simulación
    theta_ini = 0.14; %Humedad inicial (cm3/cm3)

    Dt = 1/24; %1/día Dt es por cada hora
    Dz = 15; %cm

    D_H = 24; %conversión día/hora
    A = 308; %área de irrigación m2

    Z_vect = [0:Dz:H];
    T_vect = [0:Dt:T];
    Nz = H/Dz; %número de puntos de simulación en eje z
    Nt = T/Dt; %número de puntos de simulación en eje t

    theta_v = ones(Nz+1,1)*theta_ini;
    theta_act = theta_v;
    k = 1;
    D_full = D_matrx.*(1.5*10000);
    K_full = K_matrx.*10;
    D_ini = D_full(:,1);
    K_ini = K_full(:,1);
    theta_m_pred(:,k) = theta_act;
    n_d = length(Q_0_vect);

    T_dias = zeros(length(T),1);
    for d = 1:n_d %paso del tiempo en días
    
        T_dias(d,1) = d;
        
        q0 = Q_0_vect(d,1); %m3/d Caudal de ingreso (VAR DE DECISION) se cambia el caudal en cada día (cada vez que pasa el ciclo for superior
        
        R = (q0/A)*100; % tasa de riego en (cm/d)
        
        for h = 1:D_H %contador de horas por cada día
            
            k = (d-1)*D_H + h; %contador de horas correlativas
            
            K_step = K_full(:,k);
            D_step = D_full(:,k);

            b = Pred_correct_b_PINN(theta_act,R,par_solid,dim_model,K_step,D_step);
            theta_v = theta_act;
            N = 20;
            theta_v_p_iter = zeros(Nz+1,N);
            
            for j = 1:N
                theta_v_p_iter(:,j) = Pred_correct_PINN(theta_v,b,R,par_solid,dim_model,net_k,net_D);
                theta_v = theta_v_p_iter(:,j);
            end
            theta_act = theta_v;
            theta_m_pred(:,k+1) = theta_act;
        end
    end
    theta_out = theta_m_pred;

end