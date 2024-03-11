function [out] = sigmoid_f(x)
  out = 1 ./(1+exp(-x));
end