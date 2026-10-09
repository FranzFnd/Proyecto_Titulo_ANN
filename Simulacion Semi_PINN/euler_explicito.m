function [Se_next] = euler_explicito(Se, D, K, Rj, Param)
    % Obtiene la dimensión del vector Se
    dt=Param(1);tau=Param(2);H=Param(3);dchi=Param(4);dtheta=Param(5);
    N = length(Se);
    D = (D').*10000;
    Se = extractdata(Se);
    Se = Se';
    K = (K').*10;
        % Inicializa el nuevo vector Se(:, j+1)
    Se_next = zeros(length(Se),1);
    So = (2*(Rj-K(1))*H*dchi)/(D(1)*dtheta)+Se(2);
    Se_next(1) = Se(1) + (dt * tau / H^2) *(((D(2) - D(1)) / (2 * dchi)) * ((Se(2) - So) / (2 * dchi)) + D(1) * (Se(2) - 2 * Se(1) + So) / (dchi^2) - (H / dtheta) * ((K(2) - K(1)) / (2 * dchi)));
    % Bucle sobre los nodos internos (sin incluir los extremos)
    for i = 2:N-1
        % Cálculo de los términos de la ecuación
        term1 = ((D(i+1) - D(i-1)) / (2 * dchi)) * ((Se(i+1) - Se(i-1)) / (2 * dchi));
        term2 = D(i) * (Se(i+1) - 2 * Se(i) + Se(i-1)) / (dchi^2);
        term3 = (H / dtheta) * ((K(i+1) - K(i-1)) / (2 * dchi));
        
        % Aplicación del método de Euler explícito
        Se_next(i) = Se(i) + (dt * tau / H^2) * (term1 + term2 - term3);
    end
    Se_next(N) = Se_next(N-1);

end
