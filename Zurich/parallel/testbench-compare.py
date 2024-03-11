import pandas as pd
import matplotlib.pyplot as plt

# Read the CSV file
data = pd.read_csv('compare_cvm.csv',delimiter=",")

# Sort the data based on 'size'
data = data.sort_values(by='size')

# Get unique function names
functions = data['function'].unique()
grp = data.groupby(['function', 'size'])

grp = (grp.agg({'time_us': ['mean', 'std', 'min', 'max']})
    #   .assign(time_us=lambda x:  x.time_us/1000)
      )


plt.figure(figsize=(10, 6))

for func in functions:
    func_result = grp.loc[func] 
    plt.errorbar(func_result.index.get_level_values('size'), func_result['time_us']['mean'], 
                 yerr=func_result['time_us']['std'], marker='o', label=func, 
                 linestyle='-', capsize=5, capthick=2)
    

plt.xlabel('Dimensão Matricial(NxN)')
plt.ylabel('Execution Time (ms)')
plt.title('Comparativo de Execução da Predição de Matriz de Covariancia')
plt.ticklabel_format(axis='y', style='sci',  scilimits=(0,0))
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("bench-gemm")
# plt.show()

