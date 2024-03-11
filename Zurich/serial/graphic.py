
import matplotlib.pyplot as plt 
import pandas as pd
import os
import sys
# Leitura de Dados, perfilamento
# file_path = './exec_time_omp.out'
def print_graphic(file_path):
    if not os.path.exists(file_path): 
        print("File not Found\n")
        exit(1)

    # Read File
    df = pd.read_csv(file_path, sep=',')
    
    # Create an unique call for functions
    grp = df.groupby(['function'])['time_us']

    # Create a new df with mean and std, sort from mean
    grp = grp.agg(['mean','std', 'min', 'max']).sort_values(by='mean', ascending=False)
    print(grp)
    # Time to plot
    fig, ax = plt.subplots()
    bars = ax.bar(grp.index, grp['mean'], yerr=grp['std'])

    ax.set_title('Tempo exec.: '+file_path)
    ax.set_xlabel('Etapa')
    ax.set_ylabel('Tempo Médio[us]')
    ax.legend(['Média'])

    # Rotacionar as etiquetas do eixo x para melhor legibilidade
    plt.xticks(rotation=45, ha='right')

    # Adicionar valor médio como texto acima de cada barra
    for bar in bars:
        yval = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2, yval + 0.05*yval, f'{yval:.2f}', ha='center', va='bottom')

    plt.show()

if len(sys.argv) >1:
    for i in range(1, len(sys.argv)):
        print_graphic(sys.argv[i])
else:
    print("Usage:", sys.argv[0],"<file_name>")
