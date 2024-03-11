function [D] = hs1(s,X)
  %veotriza X
  W1 = reshape(s(1:12),[3,4])' ; %W [4,3]

  v1  = W1 * X;
  y1 = sigmoid_f(v1)';
  
  D = y1;
 end