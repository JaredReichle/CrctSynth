import numpy as np

def fitcalcABCDE(sk,A,B,C,D,E):
    Nc = len(D)
    #N = len(A)
    #I = np.ones([N,1])
    
    dum = np.tile(1/(sk-A),1,Nc)
    C = C*np.transpose(dum)
    
    Yfit = C*B + D + sk*E
    
    return Yfit

