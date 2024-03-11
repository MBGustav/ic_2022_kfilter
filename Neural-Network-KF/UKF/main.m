clear all 
close all 
clc;


X = [0 0 1;
     0 1 1;
     1 0 1;
     1 1 1];

%D = [0; 1; 1; 0]; %XOR
%D = [0; 1; 1; 1]; %OR
D = [0; 0; 0; 1]; %AND 
%D = [0; 0; 1; 1]; %X



W =  randn(1,3);

%%Condicoes iniciais
N = 1000;
ns = size(W, 2);
nd = size(D, 1);
s = W';
kappa = 0;
kmax = 2*ns+1;
 
%% Modelagem do Sistema
P = 0.01*eye(ns);
R = 1.0*eye(nd);
Q = 1.0*eye(ns);
z = D;

for j=1:N
    %% Measure SigmaPoints
    [Si Ws] = SigmaPoints(s, P, kappa, kmax);
    %%  Estado de Update
     % Propagacao em fs(): Descrição do Sistema
    fSi = zeros(ns, kmax);
    for k=1:kmax
      fSi(:,k) = fs_(Si(:,k));
    end
    
    % Propagacao em hs(): Predição de Saída
    hSi = zeros(nd, kmax);
    for k=1:kmax
      for i =1:nd
        h(i,:) = hs_(Si(:,k)', X(i,:)');
      end
      hSi(:,k) = h;
    end
    
      % Transformação Uncented - UT
    [sp Pp] = UT(fSi, Ws, Q, kmax);
    [zp Pz] = UT(hSi, Ws, R, kmax);
    
    Psz = zeros(ns, nd);
    for k=1:kmax
      % Matriz de Covariancia => Psz
       Psz = Psz + Ws(k) * (fSi(:,k) - sp) * (hSi(:,k) - zp)';
    end
    
    error = z - zp;
    % Estado de Correcao
    K = Psz* Pz^(-1);
    s = sp + K*error;
    P = Pp - K*Pz*K';
    
    cost(j) = sum(abs(error));
end

%%Fase de Reconhecimento - saidas obtidas
for i=1:nd
  x = X(i,:);
  y = sigmoid_f(s'*x')
end 

plot(cost)