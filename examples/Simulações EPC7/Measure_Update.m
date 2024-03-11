function [x, Sigma] = Measure_Update(x_pre, y_pre, X_pre, Y_pre, Sigma_pre, sensor, Wc)

    % Criação de vetores e matrizes
    ids = sensor.id;
    ranges = sensor.range;
    tam = numel(ids);
    L = numel(x_pre);
    
    % Ruído do processo
    R = 0.5*eye(tam);
    
    % Preenchimento do vetor de leitura dos sensores
    for i=1:tam
        z(i) = ranges(i);
    end
  
    % Matriz de covariância entre y e y
    Pyy = zeros(tam);
    for k=1:tam
        Pyy = Pyy + Wc(k)*(Y_pre(:, k) - y_pre)*(Y_pre(:, k) - y_pre)';
    end
    Pyy = Pyy + R;
    
    % Matriz de covariância entre y e x
    Pxy = zeros(L,tam);
    for k=1:tam
        Pxy = Pxy + Wc(k)*(X_pre(:, k) - x_pre)*(Y_pre(:, k) - y_pre)';
    end
    
    % Filtragem de Kalman
    K = Pxy*inv(Pyy);
    x = x_pre + K*(z' - y_pre);
    Sigma = Sigma_pre - K*Pyy*K';
    
end