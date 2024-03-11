import numpy as np
from numpy.linalg import inv
from scipy.linalg import cholesky
import matplotlib.pyplot as plt

'''Definição da sigmoid_f'''
def sigmoid_f(x):
    return 1/(1+np.exp(-x))


'''Definicao da Transicao de Estados'''
def fs_(sp): 
    return sp


'''Definicao da Matriz de Medidas'''
def hs_(W,X):
    dot = np.dot(W, X)
    D = sigmoid_f(dot)
    return D

def hs1(s, X):  
    W1 = s[ 0:12].reshape(4, 3)
    W2 = s[12:16].reshape(1, 4)

    v1 = np.matmul(W1, X)
    y1 = sigmoid_f(v1)
    return y1

def hs2(s, X):
    W1 = s[ 0:12].reshape(4, 3)
    W2 = s[12:16].reshape(1, 4)
    '''Propagacao primeira camada'''
    v1 = np.matmul(W1, X)
    y1 = sigmoid_f(v1)
    '''Propaga para segunda camada'''
    v = np.matmul(W2, y1)
    y = sigmoid_f(v)
    return y




'''Definicao de Funcoes do Filtro do Kalman'''
def SigmaPoints(xm, P, kappa, kmax):
    n = xm.shape[0] # Considerando vetor-coluna!
    Xi = np.zeros((n,kmax))
    Ws = np.zeros((kmax,1))

    # Decomp. cholesky: U x U.T = (n+kp)*P
    U = cholesky((n+kappa)*P)
    Ws[0] = kappa / (n+kappa)
    Xi[:,0] = xm.flatten()
    for k in range(n):
        Xi[:, k+1] = xm.flatten() + U[k, :].T
        Ws[k+1] = 1 / (2*(n+kappa))

    for k in range(n): 
        Xi[:,n+k+1] = xm.flatten() - U[k, :].T
        Ws[n+k+1] = 1 / (2*(n+kappa))

    return Xi, Ws

def UT(Xi, W, noiseCov): 
    n, kmax = Xi.shape
    xm = np.zeros((n,1))
    xcov = np.zeros((n,n))

    '''Acumulacao do vetor xm'''
    for k in range(kmax):
        xm += W[k] * Xi[:, k].reshape((n,1))

    '''Acumulacao da matriz da Cov. de Medidas'''
    for k in range(kmax): 
        vec = (Xi[:,k].reshape((n,1)) - xm)
        xcov += W[k] * np.matmul(vec , vec.T)  
    xcov += noiseCov
    return xm, xcov

'''Declaracao de Entrada e Saida Esperada'''
X = np.array([[0, 0, 1],
              [0, 1, 1],
              [1, 0, 1],
              [1, 1, 1]])

# D = np.array([[0, 0, 0, 1]]).T #AND
D = np.array([[0, 1, 1, 1]]).T #OR
# D = np.array([[0, 0, 1, 1]]).T #X 

'''XOR - nao é possivel treinar com uma camada'''
# D = np.array([[0, 1, 1, 0]]).T #XOR


'''Condicoes Iniciais'''
N = 300
nd = D.shape[0]
ni = 3
nh = 4 
no = 1
NLayers = 2
cost = np.ones((N,))
z = D
t = np.arange(N)


W1 = 2*np.random.random((nh,ni)) - 1 
W2 = 2*np.random.random((no,nh)) - 1 



ns1 = W1.size
ns2 = W2.size
ns = ns1+ns2

s = np.concatenate((np.hstack(W1), np.hstack(W2)), axis=None)


'''Modelagem do Sistema'''
kappa = 0
kmax = 2*ns+1
P = 0.01* np.eye(ns)
R2 = 100   * np.eye(nd*no)
R1 = 0.1   * np.eye(nd*nh)
Q = 0.0001 * np.eye(ns)


for j in range(N):    
    
    '''Processo de Sigma Points'''
    Si, Ws = SigmaPoints(s, P, kappa, kmax)
    fSi = np.zeros((ns,kmax))

    '''Estado de Update'''
    '''1. Propagacao em f(s)'''
    for k in range(kmax):
        fSi[:, k] = fs_(Si[:, k])
    

    for layer in range(NLayers):
        zaux1 = np.zeros((nh, nh))
        if layer == 0:
            nn = nh
            R = R1
        else:
            nn = no
            R = R2

        '''2. Propagacao em h(s)'''
        hSi = np.zeros((nd*nn, kmax))
        h = np.zeros((nd,nn))

        for k in range(kmax):
            for i in range(nd): 
                if layer == 0:
                    h[i, :] = hs1(Si[:, k], X[i, :])
                else:
                    h[i, :] = hs2(Si[:, k], X[i, :])
        
            if layer == 0: 
                hm, hn = h.shape
                hSi[:, k] = h.flatten()
            else :
                hSi[:, k] = h.flatten()
        
        '''Transformacao Unscented - UT'''
        sp, Pp = UT(fSi, Ws, Q)
        zp, Pz = UT(hSi, Ws, R)

        Psz = np.zeros((ns, nd*nn))
        for k in range(kmax):
            v1 = (fSi[:, k].reshape((ns, 1)) - sp).reshape((ns,1))
            v2 = (hSi[:, k].reshape((nd*nn, 1)) - zp).reshape((nd*nn,1))
            Psz += Ws[k] * np.matmul(v1, v2.T)
        if layer == 0: 
            for k in range(nd):
                zaux1[:, k] = np.matmul(W2.T, D[k])
            z = zaux1.T.reshape([nh*nh, 1])
        else:
            z = D

        '''Estado de Correcao'''
        error = z - zp
        inv_Pz = inv(Pz)
        K = np.matmul(Psz ,inv_Pz)
        s = sp + np.matmul(K, error)
        P = Pp - np.matmul(K, np.matmul(Pz,K.T))
        W1 = s[:12].reshape(4, 3)
        W2 = s[12:].reshape(1, 4)

    cost[j] = np.sum(np.abs(error))
    print(cost[j])


y = np.zeros(nd)
for i in range(nd): 
    x = X[i,:]
    W1 = s[ 0:12].reshape(4, 3)
    W2 = s[12:16].reshape(1, 4)

    v1 = np.matmul(W1,x)
    y1 = sigmoid_f(v1)
    v = np.matmul(W2, y1)
    y[i] = sigmoid_f(v)




plt.figure()
plt.plot(t, cost)
plt.show()
