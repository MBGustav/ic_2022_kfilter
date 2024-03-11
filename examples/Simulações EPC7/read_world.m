function [landmarks] =  read_world(world)
    % Lê a definição do mundo e retorna uma lista de pontos de referência, 
    % ou seja o 'mapa'
    %
    % O dicionário retornado contém uma lista de pontos de referência, cada
    % um com a seguinte informação: [id, x, y]

    landmarks(:,1) = world.VarName1(1:end);
    landmarks(:,2) = world.VarName2(1:end);
    landmarks(:,3) = world.VarName3(1:end);
end