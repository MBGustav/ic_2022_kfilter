clear all
syms w1 w2 w3
syms x1 x2 x3

W = [w1;w2;w3];
X = [x1 x2 x3];


D = tanh(X*W);

Hl(1,1) = diff(D,w1);
Hl(1,2) = diff(D,w2);
Hl(1,3) = diff(D,w3);

Hl
