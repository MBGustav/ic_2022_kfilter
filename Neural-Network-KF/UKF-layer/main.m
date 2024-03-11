clear all 
close all 
clc;


X = [0 0 1;
     0 1 1;
     1 0 1;
     1 1 1];

D = [0; 1; 1; 0]; %XOR
%D = [0; 1; 1; 1]; %OR
%D = [0; 0; 0; 1]; %AND 
%D = [0; 0; 1; 1]; %X



%% Gera W random entre -1 e 1
W1 = 2*rand(4,3)-1
W2 = 2*rand(1,4)-1

%%Condicoes iniciais
N = 5000;

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

kappa = 0;
kmax = 2*ns+1;
 
%% Modelagem do Sistema
P = 0.01*eye(ns);
%R = 1.0*eye(nd);
R2 = 1.0*eye(nd*no);
R1 = 100*eye(nd*nh);
Q = 0.001*eye(ns);


for j=1:N
    %% Measure SigmaPoints
    [Si Ws] = SigmaPoints(s, P, kappa, kmax);
    %%  Estado de Update
     % Propagacao em fs(): Descrição do Sistema
    fSi = zeros(ns, kmax);
    for k=1:kmax
      fSi(:,k) = fs_(Si(:,k));
    end
    
    
    Nlayer = 2;
    for layer = 1: Nlayer
        if layer == 1
            nn = nh; % número de neurônios
            R = R1; 
        else 
            nn = no;
            R = R2;
        end
        % Propagacao em hs(): Predição de Saída
        hSi = zeros(nd*nn, kmax);
        h = zeros(1,nn);
        for k=1:kmax
          for i =1:nd
            if layer == 1
                h(i,:) = hs1(Si(:,k)', X(i,:)');
            else
                h(i,:) = hs2(Si(:,k)', X(i,:)');
            end
          end
          if layer == 1
              [hm,hn]=size(h);
              hSi(:,k) = reshape(h',[1,hm*hn]);
          else
            hSi(:,k) = h;
          end
        end

      
      % Transformação Uncented - UT
        [sp Pp] = UT(fSi, Ws, Q, kmax);
        [zp Pz] = UT(hSi, Ws, R, kmax);

        Psz = zeros(ns, nd*nn);
        for k=1:kmax
          % Matriz de Covariancia => Psz
           Psz = Psz + Ws(k) * (fSi(:,k) - sp) * (hSi(:,k) - zp)';
        end
        if layer == 1
            for k= 1:nd
                zaux1(:,k) = W2*D(k);
            end
            z = reshape(zaux1,[16,1]);
        else
            z = D;
        end

        error = z - zp;
        % Estado de Correcao
        K = Psz* Pz^(-1);
        s = sp + K*error;
        P = Pp - K*Pz*K';
        W1 = reshape(s(1:12),[3,4])' ; %W [4,3]
        W2 = reshape(s(13:16),[4,1])'; %w [1,4]
    end
    cost(j) = sum(abs(error));
    disp(cost(j));
    
%     figure(1)
%     plot(cost)
end

%%Fase de Reconhecimento - saidas obtidas
for i=1:nd
  x = X(i,:)';
  %y = sigmoid_f(s'*x')
  
  W1 = reshape(s(1:12),[3,4])' ; %W [4,3]
  W2 = reshape(s(13:16),[4,1])'; %w [1,4]
  
  v1  = W1 * x;
  y1 = sigmoid_f(v1);
  v = W2*y1;
  y(i) = sigmoid_f(v);

end
  
  %%Peso
  s;
  

plot(cost)
