function [out] = Hs(W,X)
%Sem bias!
v = W*X;
y = sigmoid_f(v);
aux  = y*(1- y)';

out(1,1) = aux*X(1);
out(1,2) = aux*X(2);
out(1,3) = aux*X(3);
end