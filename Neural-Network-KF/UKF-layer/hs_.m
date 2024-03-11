function [D] = hs_(s,X)
  %veotriza X
  W1 = reshape(s(1:12),[3,4])' ; %W [4,3]
  W2 = reshape(s(13:16),[4,1])'; %w [1,4]
  
  v1  = W1 * X;
  y1 = sigmoid_f(v1)';
  
  v = W2*y1;
  y = sigmoid_f(v);
 
  D = y;
 end