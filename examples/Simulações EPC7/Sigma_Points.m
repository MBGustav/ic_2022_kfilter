function [Xi,Wm,Wc] = Sigma_Points(xm,Sigma)

    alpha  = 1e-2;
    L      = numel(xm);
    k      = 0;
    beta   = 1;
    lambda = alpha^2*(L + k) - L;

    % Criação das matrizes
    Xi = zeros(L,2*L+1);
    Wm = zeros(L,1);
    Wc = zeros(L,1);

    % Vetor Sigma0 e pesos relacionados
    Xi(:,1) = xm;
    Wm(1)   = lambda/(L + lambda);
    Wc(1)   = lambda/(L + lambda) + (1 - alpha^2 + beta);

    % Decomposição de Cholesky % U'*U = (L+lambda)*P
    U = chol((L+lambda)*Sigma); 

    % Preenchimento da matriz de Sigma Points e pesos
    for k=1:L
        Xi(:,k+1) = xm + U(:, k);
        Wm(k+1)   = 1 / (2*(L+lambda));
        Wc(k+1)   = 1 / (2*(L+lambda));
    end

    for k=1:L
        Xi(:,L+1+k) = xm - U(:, k);
        Wm(L+k+1)   = 1 / (2*(L+lambda));
        Wc(L+k+1)   = 1 / (2*(L+lambda));
    end
    
end