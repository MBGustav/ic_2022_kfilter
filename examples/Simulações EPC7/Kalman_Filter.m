% Limpeza de variáveis, janelas e área de trabalho
clear all
close all
clc

load sensor_data.mat
load world.mat

disp('Reading landmark positions')
landmarks = read_world(world);

disp('Reading sensor data')
sensor_readings = read_sensor_data(sensordata);

% Valores iniciais da distribuição normal
x = [0, 0, 0]';
Sigma = eye(3);
     
 % Ruído do processo
 Q = [0.2, 0.0, 0.0;
      0.0, 0.2, 0.0;
      0.0, 0.0, 0.02];

% Intervalo da janela de animação
map_limits = [-3, 13, -2, 12];

% Tamanho de estados e vetor de dados
tam = numel(x);
N = length(sensor_readings);

% Variáveis para salvar
muv_pre = zeros(tam,N);
muv = zeros(tam,N);
timestepv = zeros(tam,N);
    
for timestep = 1:N
    
%     % Traça o estado atual
%     Plot_State(x, Sigma, landmarks, map_limits);
    
    % Cálculo dos Sigma Points
    [Xi,Wm,Wc] = Sigma_Points(x, Sigma);

    % Informações de sensores e odometria
    sensor = sensor_readings(timestep).sensor;
    odometry = sensor_readings(timestep).odometry;
    
    % Time Update
    [x_pre, y_pre, X_pre, Y_pre, Sigma_pre] = Time_Update(sensor, odometry, landmarks, x, Xi, Wm, Wc, Q);

    % Measurerement Update Equations
    [x, Sigma] = Measure_Update(x_pre, y_pre, X_pre, Y_pre, Sigma_pre, sensor, Wc);
         
    muv_pre(:,timestep) = x_pre;
    muv(:,timestep) = x;
    timestepv(timestep) = timestep;
       
end
  

%% Gráficos

    figure(1)
    plot(muv(1,:),muv(2,:),'r','LineWidth',1.5)
    title('\it Trajetória do Robô no Plano','FontSize',14)
    xlim([-1 11])
    ylim([-1 11])
    grid on
    grid minor
    xlabel('\textit{Deslocamento em x}','FontName','Times New Roman','FontSize',14,'Interpreter','latex')
    ylabel('\textit{Deslocamento em y}','FontName','Times New Roman','FontSize',15,'Rotation',90,'Interpreter','latex')
    set(gca,'XTick',(-1:11),'YTick',(-1:11),'Box','on','XMinorTick','on','YMinorTick','on')
    set(gcf,'Position',[300 100 700 500],'PaperPositionMode','auto')
    
    figure(2)
    plot(muv(3,:)*180/pi,'r','LineWidth',1.5)
    title('\it Orientação do Robô no Plano','FontSize',14)
    xlim([0 330])
    ylim([-300 150])
    grid on
    grid minor
    xlabel('\textit{Steps} [n]','FontName','Times New Roman','FontSize',14,'Interpreter','latex')
    ylabel('\textit{Angulo} ($^{\circ}$)','FontName','Times New Roman','FontSize',15,'Rotation',90,'Interpreter','latex')
    set(gca,'XTick',(0:30:330),'YTick',(-300:50:150),'Box','on')
    set(gca,'XTickLabel',{'0','30','60','90','120','150','180','210','240','270','300','330'},'XMinorTick','on');
    set(gca,'YTickLabel',{'-300','-250','-200','-150','-100','-50','0','50','100','150'},'YMinorTick','on');
    set(gcf,'Position',[300 100 700 500],'PaperPositionMode','auto')