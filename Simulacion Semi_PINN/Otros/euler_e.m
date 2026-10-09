function [Se_next] = euler_e(Se, D, K, Rj, Param)
    % Obtiene la dimensión del vector Se
    dt=Param(1);tau=Param(2);H=Param(3);dchi=Param(4);dtheta=Param(5);
    % D = (D).*1000;
    % Se = extractdata(Se);
    % Se = Se';
    % K = (K).*100;
    N = length(Se);
        % Inicializa el nuevo vector Se(:, j+1)
    Se_next = zeros(length(Se),1);
    K11 = K(1,1);D11 = D(1,1);Se21 = Se(2,1);Se11 = Se(1,1);D21=D(2,1);K21 = K(2,1);
    So = (2*(Rj-K11)*H*dchi)/(D11*dtheta)+Se21;
    Se_next(1,1) = Se11 + (dt * tau / H^2) *(((D21 - D11) / (2 * dchi)) * ((Se21 - So) / (2 * dchi)) + D11 * (Se21 - 2 * Se11 + So) / (dchi^2) - (H / dtheta) * ((K21 - K11) / (2 * dchi)));
    % Bucle sobre los nodos internos (sin incluir los extremos)
    for i = 2:N-1
        % Cálculo de los términos de la ecuación
        term1 = ((D(i+1,1) - D(i-1,1)) / (2 * dchi)) * ((Se(i+1,1) - Se(i-1,1)) / (2 * dchi));
        term2 = D(i,1) * (Se(i+1,1) - 2 * Se(i,1) + Se(i-1,1)) / (dchi^2);
        term3 = (H / dtheta) * ((K(i+1,1) - K(i-1,1)) / (2 * dchi));
        
        % Aplicación del método de Euler explícito
        Se_next(i,1) = Se(i,1) + (dt * tau / H^2) * (term1 + term2 - term3);
    end
    Se_next(N) = Se_next(N-1);

end
