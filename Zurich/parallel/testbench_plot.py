import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Leitura dos dados do arquivo CSV
file_path = 'testbench_CovMatrix.csv'  # Substitua 'seu_arquivo.csv' pelo caminho real do seu arquivo
data = pd.read_csv(file_path)

# Filtro dos dados com base no tamanho (size = 15)
filtered_data = data[data['size'] == 15]

# Exclusão dos dois valores máximos de cada grupo
filtered_grouped_data = filtered_data.groupby(['type_execution', 'padding', 'device_type']).apply(lambda x: x[x['time_us'] < x['time_us'].nlargest(2).min()])

# Redefinir o índice para evitar ambiguidade
filtered_grouped_data = filtered_grouped_data.reset_index(drop=True)

# Cálculo da média dos tempos de execução
average_execution_times = filtered_grouped_data.groupby(['type_execution', 'padding', 'device_type']).agg({'time_us': 'mean'}).unstack(level=['padding', 'device_type']).fillna(0)

# Verificação e preenchimento de grupos vazios
if average_execution_times.empty:
    print("Não há dados suficientes para gerar o gráfico.")
else:
    # Preencher níveis ausentes no índice após o unstack e substituir NaN por 0
    reshaped_data = average_execution_times

    # Criação do gráfico de barras
    fig, ax = plt.subplots()

    bar_width = 0.2
    bar_positions = np.arange(len(reshaped_data.index))

    # Barras para o padding 'Y' no device 'CPU'
    if ('time_us', 'Y', 'CPU') in reshaped_data.columns:
        bar_y_cpu = ax.bar(bar_positions - bar_width, reshaped_data[('time_us', 'Y', 'CPU')], bar_width, label='Padding Y - CPU')

    # Barras para o padding 'Y' no device 'GPU'
    if ('time_us', 'Y', 'GPU') in reshaped_data.columns:
        bar_y_gpu = ax.bar(bar_positions, reshaped_data[('time_us', 'Y', 'GPU')], bar_width, label='Padding Y - GPU')

    # Barras para o padding 'N' no device 'CPU'
    if ('time_us', 'N', 'CPU') in reshaped_data.columns:
        bar_n_cpu = ax.bar(bar_positions + bar_width, reshaped_data[('time_us', 'N', 'CPU')], bar_width, label='Padding N - CPU')

    # Barras para o padding 'N' no device 'GPU'
    if ('time_us', 'N', 'GPU') in reshaped_data.columns:
        bar_n_gpu = ax.bar(bar_positions + 2 * bar_width, reshaped_data[('time_us', 'N', 'GPU')], bar_width, label='Padding N - GPU')

    # Configurações do gráfico
    ax.set_xlabel('Tipo de Execução')
    ax.set_ylabel('Tempo de Execução (us)')
    ax.set_title('Tempo de Execução por Tipo, Padding e Device (Size = 15)')
    ax.set_xticks(bar_positions + bar_width/2)
    ax.set_xticklabels(reshaped_data.index)
    ax.legend()

    # Exibir o gráfico
    plt.show()
