import numpy as np 
import matplotlib.pyplot as plt
from numpy.linalg import inv
from scipy.linalg import cholesky

def fs_(s,X):
    n = 3
    A = np.eye(3)
    return np.matmul(A,s)

def hs_(s,X):
    W = s
    dot = np.dot(X,W)
    D = np.tanh(dot)
    return D

def Hs(s,X):
    out = np.zeros((3,))
    '''Produto Interno entre pesos e Entradas'''
    # dot = np.dot(s,X)
    # sum = 1-np.tanh(dot)**2
    sum = 1 - np.tanh(s[0]*X[0] + s[1]*X[1]+s[2]*X[2])**2
    
    out[0] = X[0] * sum
    out[1] = X[1] * sum
    out[2] = X[2] * sum
    return out

def Hjacob(s,X):
    Hjacobian = np.zeros((3,1))
    # dot = np.dot(X,s)
    dot_prod = s[0]*X[0] + s[1]*X[1]+s[2]*X[2]
    sum =1 - np.tanh(dot_prod)**2

    Hjacobian[0] = s[0] * sum
    Hjacobian[1] = s[1] * sum
    Hjacobian[2] = s[2] * sum
    
    return Hjacobian

def SigmaPoints(xm, P, kp): 
    n = xm.size
    Xi = np.zeros((n,2*n+1))
    W = np.zeros((2*n+1,))

    Xi[:,0] = xm.flatten()
    W[0] = kp / (n+kp)
    
    # Scalar multiplication
    U = cholesky((n+kp) * P)

    for k in range(n):
        Xi[:, k+1] =  xm.flatten() + U[:,k]
        W [k+1] = 1/(2*(n+kp))

    for k in range(n):
        Xi[:,n+k+1] = xm.flatten() - U[:,k]
        W [n+k+1] = 1/(2*(n+kp))
    
    return Xi, W



D = np.array([[0, 0, 1, 1]]).T
X = np.array([[0, 0, 1],
              [0, 1, 1],
              [1, 0, 1],
              [1, 1, 1]])

# W = np.random.random((3,1)) -1 
W = np.array([[-0.47399, 0.73467, 1.5124]]).T

'''Initial Conditions'''
N = 1
cost = np.zeros((N,))
t = np.arange(N)


ns = W.shape[0]
nd = D.shape[0]
s = W
P = np.eye(ns)
R = 5*np.eye(nd)
Q = 0.01*np.eye(ns)

H = np.zeros((nd, ns))
h = np.zeros((nd, 1))

for i in range(N):
    
    '''SigmaPoints'''
    Xi, sW = SigmaPoints(s, P, 0)
    kmax = 2*ns + 1
    fXi = np.zeros((ns, kmax))
    hXi = np.zeros((nd, kmax)) #??

    '''Predict'''
    for k in range(kmax):
        fXi[:, k] = Hs(s, Xi[:, k])
        hXi[:, k] = hs_(s, Xi[:, k])
        # print("Hs", Hs(s, Xi[:, k]))
        # print("hs_",hs_(s, Xi[:, k]))
    print("fXi",fXi)
    print("hXi",hXi)

