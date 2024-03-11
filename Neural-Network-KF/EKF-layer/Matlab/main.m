clear all
close all
clc

X = [0 0 1;
     0 1 1;
     1 0 1;
     1 1 1];
 
D = [0; 1; 1; 1];
%% Gera W random entre -1 e 1
W1 = 2*rand(4,3)-1
W2 = 2*rand(1,4)-1

%%Condicoes iniciais
N = 10;
NLayer = 2
ni = 4; %  Número de entradas
nh = 4; %  Número de neurônios na camada escondida
no = 1; %  Número de neurônios na camada de saída

ns1 = size(W1,1)*size(W1,2);
ns2 = size(W2,1)*size(W2,2);
ns = ns1 + ns2; 
nd = size(D,1);
%Ordenamos s por ordem de neuronio s = [W1,1..4 ; W2.,1..4 ; ... ] 
s = [reshape(W1',[ns1, 1]); reshape(W2',[ns2, 1])]
%W1 = reshape(s(1:12),[3,4])'
%W2 = reshape(s(13:16),[4,1])'
 
%% Modelagem do Sistema
P = 0.01*eye(ns);
%R = 1.0*eye(nd);
R2 = 1.0*eye(nd*no);
R1 = 0.1*eye(nd*nh);
Q = 0.00001*eye(ns);


for k = 1:N
    %%Predict
    sp = s;
    Pp = P + Q;
    %%Update
    z = D;
    for layer=1:Nlayer
      if layer == 1
        nn = nh;
        R = R1;
      else
        nn = no;
        R = R2; 
      endif
      
      h = zeros((1,nn));
      for i =1:nd
        if layer == 1
          H(i,:) = Hs1(s',X(i,:)');
          h(i,:) = hs_(s',X(i,:)'); % primeira camada
        else
          H(i,:) = Hs2(s',X(i,:)');
          h(i,:) = hs2_(s',X(i,:)'); % segunda camada
        endif
      end  
      
    endfor
    
      
    error = z - h;
    A = H*P*H' + R;
    K  = Pp*H'*A^(-1);
    s = sp +K*error;
    P = (eye(ns) -K*H)*Pp;
    cost(k) = sum(abs(error));
end

plot(cost)

 
