import numpy as np
import numpy.linalg as lin

def pr2ss(SER):
    R = SER.R
    poles = SER.poles
    Nc = len(SER.D)
    N = len(R[0,0,:])
    C = np.zeros([Nc,Nc*N])
    A = np.zeros([Nc*N,1], dtype = 'complex128')
    B = np.zeros([Nc*N,Nc], dtype = 'complex128')
    #poles = lin.diag(lin.diag(poles))
    for m in range(0,N):
        Rdum = np.squeeze(R[:,:,m])
        for n in range(0,Nc):
            ind = (n-1)*N+m
            C[:,ind] = Rdum[:,n]
    
    for n in range(0,Nc):
        A[(n-1)*N+1:n*N] = poles
        B[(n-1)*N+1:n*N,n] = np.ones([N,1])
    A = lin.diag(A)
    
    SER.A = A
    SER.B = B
    SER.C = C
    
    return SER

