function [theta_v_p] = Pred_correct(theta_v,b,R,par_solid, par_model,dim_model)

%PRED_CORRECT Summary of this function goes here
%   Detailed explanation goes here
theta_s = par_solid(1); theta_r = par_solid(2); Ks = par_solid(3); 

alpha = par_model(1); n = par_model(2); method = par_model(3); 
Nz = dim_model(1); Dz = dim_model(2); Dt = dim_model(3); 

%cálculos para A y c
S = S_theta(theta_v,theta_r,theta_s,Nz);
h = h_theta(S,n,alpha,Nz,method);
K = K_theta(S,Ks,n,Nz,method);
dtheta_dh = dtheta_dh_theta(h,theta_s,theta_r,n,alpha,Nz,method);
D = D_theta(K,dtheta_dh,Nz);
D0 = D(1);
K0 = K(1);
theta_0 = theta_v(1);

[theta_m1, D_m1, K_m1] = theta_m1_theta(theta_0,D0,K0,Dz,R,n,alpha,Ks,theta_r,theta_s,method);




%Armar las matrices y vectores asociados 
w1 = Dt/(2*Dz^2);
w2 = Dt/(4*Dz);


%MATRIZ A

A = zeros(Nz+1);
%b = zeros(Nz+1,1);
c = zeros(Nz+1,1);



%nodo superior
i = 1; 
Dff = (D(i+1) + D(i))/2;
Dbf = (D(i) + D_m1)/2;

%Df = (D(i+1) + D(i))/2;
%Db = (D(i) + D_m1)/2;

A(i,i) = 1+w1*(Dff + Dbf);
A(i,i+1) = -w1*Dff;

%b(i,1) = theta_v(i) + w1*(Df*(theta_v(i+1)-theta_v(i)) - Db*(theta_v(i)-theta_m1)) - w2*(K(i+1) - K_m1);
c(i,1) = -w2*(K(i+1) - K_m1) + w1*Dbf*theta_m1;

%nodos interiores
for i = 2:Nz
    Dff = (D(i+1) + D(i))/2;
    Dbf = (D(i) + D(i-1))/2;
    
    %Df = (D(i+1) + D(i))/2;
    %Db = (D(i) + D(i-1))/2;
    
    A(i,i) = 1+w1*(Dff + Dbf);
    A(i,i+1) = -w1*Dff;
    A(i,i-1) = -w1*Dbf;
    
%    b(i,1) = theta_v(i) + w1*(Df*(theta_v(i+1)-theta_v(i)) - Db*(theta_v(i)-theta_v(i-1))) - w2*(K(i+1) - K(i-1));
    c(i,1) = -w2*(K(i+1) - K(i-1));
    
end

%nodo inferior
i = Nz+1;
Dff = D(i);
Dbf = (D(i) + D(i-1))/2;
%Df = D(i);
%Db = (D(i) + D(i-1))/2;

A(i,i) = 1+w1*(Dff + Dbf);
A(i,i-1) = -w1*Dbf;

%b(i,1) = theta_v(i) + w1*(Df*(0) - Db*(theta_v(i)-theta_v(i-1))) - w2*(K(i) - K(i-1)); %ojo --> 0 = theta_nz+1 - theta_Nz y K_Nz+1 = K_Nz por CB en el nodo inferior

c(i,1) = -w2*(K(i) - K(i-1)) + w1*Dff*theta_v(i);

B = c+b;
%intento 1
theta_v_p = A\B; %OJO  --> mejorar la inversión tridiagonal con diagonal principal dominante --> Thomas (desc LU)

%intento 2 --> invertir una matriz, descomponiendola con factorización LU
%[L,U,P] = lu(A);
%theta_v_p = U\(L\(P*B));

for i=1:Nz+1
    if theta_v_p(i,1) > theta_s
        theta_v_p(i,1) = theta_s;
    end
    
    if theta_v_p(i,1) < theta_r
        theta_v_p(i,1) = theta_r;
    end
    
end

