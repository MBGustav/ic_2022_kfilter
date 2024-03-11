function [Xi Ws] = SigmaPoints(xm, P, kappa, kmax);
  
  n = numel(xm);
  Xi = zeros(n, kmax);
  Ws = zeros(kmax, 1);
  
  U = chol((n+kappa)*P);
  Ws(1) = kappa / (n+kappa);
  Xi(:,1) = xm;
  
   for k=1:n
    Xi(:, k+1) = xm + U(k, :)';
    Ws(k+1)     = 1 / (2*(n+kappa));
  end

  for k=1:n
    Xi(:, n+k+1) = xm - U(k, :)';
    Ws(n+k+1)      = 1 / (2*(n+kappa));
  end

end