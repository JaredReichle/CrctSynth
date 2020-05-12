import numpy as np
import numpy.matlib as mat

def fitcalcABCDE(sk,A,B,C,D,E):
    Nc = len(D)
    #N = len(A)
    #I = np.ones([N,1])
    
    dum1 = (sk-A)
    dum2 = 1/dum1
    dum = mat.repmat(dum2,Nc,1)
    C = C*dum
    
    tmp1 = np.matmul(C,B)
    tmp2 = sk*E
    tmp3 = tmp1 + D + tmp2
    Yfit = tmp3
    #Yfit = C*B + D + sk*E
    
    return Yfit

