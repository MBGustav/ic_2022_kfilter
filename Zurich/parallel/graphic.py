
import matplotlib.pyplot as plt 
import pandas as pd
import os
import sys
# Leitura de Dados, perfilamento
def print_graphic(file_path):
    if not os.path.exists(file_path): 
        print("File not Found\n")
        exit(1)

    # Read File
    df = pd.read_csv(file_path, sep=',')
    
    # This remove max and min from each function....
    df = df.groupby('function').apply(lambda x: x[(x['time_us'] != x['time_us'].max()) & (x['time_us'] != x['time_us'].min())]).reset_index(drop=True)
    # Doing a second time to remove again ?
    df = df.groupby('function').apply(lambda x: x[(x['time_us'] != x['time_us'].max()) & (x['time_us'] != x['time_us'].min())]).reset_index(drop=True)
    
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

    # Rotate 45 to show in plotting
    plt.xticks(rotation=45, ha='right')

    #Include medium val for each bar
    for bar in bars:
        yval = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2, yval + 0.05*yval, f'{yval:.2f}', ha='center', va='bottom')
    plt.tight_layout()
    plt.show()

if len(sys.argv) >1:
    for i in range(1, len(sys.argv)):
        print_graphic(sys.argv[i])
else:
    print("Usage:", sys.argv[0],"<file_name>")
    print("In this file, the script makes a comparative with only time execution.")
    print("For more detailed benchmark, with a third parameter (size, for example), use: ./testbench-compare")
