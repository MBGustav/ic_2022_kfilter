function [D] = hs_(W,X)
v = W*X;
    D = sigmoid_f(v);
end