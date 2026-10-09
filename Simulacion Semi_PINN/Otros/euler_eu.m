function [Se_next] = euler_eu(Se, D, K, Rj, Param)
    
    dt=Param(1);tau=Param(2);H=Param(3);dchi=Param(4);dtheta=Param(5);
    N = length(Se);
    Se = extractdata(Se);
    Se = Se';
    D = (D').*1000;
    K = (K').*10;
   
    % Inicializa el nuevo vector Se(:, j+1)
    Se_next = zeros(length(Se),1);
    %So = (2*H*dchi)*(K(1,1)-Rj)/(dtheta*D(1,1)); % CI Centrada
    So = Se(2,1) + (H*dchi)*(-K(1,1)+Rj)/(dtheta*D(1,1)); % CI adelanto
    
    % Situacion inicial
    Se_next(1,1) = Se11 + (dt * tau / H^2) *(((D(i+1) - D(i-1)) / (2 * dchi)) * ((Se(i+1) - So) / (2 * dchi)) + D(i) * (Se(i+1) - 2 * Se(i) + So) / (dchi^2) - (H / dtheta) * ((K(i+1) - K(i-1)) / (2 * dchi)));
    % Bucle sobre los nodos internos (sin incluir los extremos)
    for i = 2:N-1
        % Cálculo de los términos de la ecuación
        term1 = ((D(i+1,1) - D(i-1,1)) / (2 * dchi)) * ((Se(i+1,1) - Se(i-1,1)) / (2 * dchi));
        term2 = D(i,1) * (Se(i+1,1) - 2 * Se(i,1) + Se(i-1,1)) / (dchi^2);
        term3 = (H / dtheta) * ((K(i+1,1) - K(i-1,1)) / (2 * dchi));
        
        % Aplicación del método de Euler explícito
        valant = Se(i,1);
        Se_next(i,1) = valant + (dt * tau / H^2) * (term1 + term2 - term3);
    end
    Se_next(N,1) = Se_next(N-1,1);

end
