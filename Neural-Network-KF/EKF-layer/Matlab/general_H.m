function [H] = general_H(U, X, list_W)
  % U       - valores de entrada da rede neural
  % X       - valores de saida esperada da rede neural
  % list_W  - lista de matriz de pesos
  % k_total: total de camadas a serem dessenvolvidas
  % 
  
  %calcula H1:
  #W1 = list_W{1};
  #v1 = W1 * U
  #H = {sigmoid_f(v1)};
  nd = size(X,1);
  k_total = size(list_W, 2);
  for k=1:k_total
    Wk = list_W{k};
    if k==1  
      y = sigmoid_f(Wk * U);
    else
      y = sigmoid_f(Wk * prev_y);
    end
    for i
    H{k} = y;
    
    prev_y = y;
  endfor
  
end
