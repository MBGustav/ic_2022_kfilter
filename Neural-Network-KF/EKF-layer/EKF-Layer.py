import numpy as np
from scipy.linalg import block_diag

def jacobian():





    pass

'''funcao sigmoid'''
def sigmoid(x): 
    return 1/(1+np.exp(-x))
'''derivada da sigmoid'''
def dx_sigmoid(x):
    return sigmoid(x)*(1-sigmoid(x))

l = np.array([])
y = np.array([])
u = np.array([])


nl = l.size(l)
ny = y.size(y) 
nu = u.size()
W = [2*(2*np.random.sample((nl, nu+1))-1),
     2*(2*np.random.sample((ny, nl+1))-1)]


print()

