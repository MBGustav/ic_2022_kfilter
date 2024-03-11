function [D] = hs_(s,X)
  v = s * X;
  D = sigmoid_f(v);
 end