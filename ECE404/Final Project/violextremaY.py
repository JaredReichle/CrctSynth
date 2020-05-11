import numpy as np
import numpy.linalg as lin
from math import pi
from fitcalcABCDE import fitcalcABCDE
from intecheig import intercheig
from rot import rot

def violextremaY(SERflag,wintervals,A,B,C,D,colinterch):
    s_pass = []
    g_pass = []
    #ss = []
    if len(wintervals) == 0:
        s_pass = []
        return
    
    if SERflag == 1:
        SERA = A
        SERC = C
        SERD = D
        Nc = len(SERD)
        N = len(SERA)
        A = np.zeros([Nc*N,1])
        for col in range(0,Nc):
            A[(col-1)*N+1:col*N] = SERA
        B = np.zeros([Nc*N,Nc])
        for col in range(0,Nc):
            B[(col-1)*N+1:col*N,col] = np.ones([N,1])
        C = np.zeros([Nc*N,Nc])
        for row in range(0,Nc):
            for col in range(0,Nc):
                C[row,(col-1)*N+1:col*N] = np.transpose(np.squeeze(SERC[row,col,:N]))
        D = SERD
        A = lin.diag(A)
        
    s = []
    EE = []
    Nc = len(D)
    g_pass = 1e16
    smin = 0
    for m in range(0,len(wintervals[:,0])):
        Nint = 21
        w1 = wintervals[m,0]
        if wintervals[m,1] == 1e16: #violation extends to inf freq
            w2 = 2*pi*1e16
        else:
            w2 = wintervals[m,1]
        s_pass1 = 1j*np.linspace(w1,w2,Nint)
        if w1 == 0:
            s_pass2 = 1j*np.logspace(-8,np.log10(w2),Nint)
        else:
            s_pass2 = 1j*np.logspace(np.log10(w1),np.log10(w2),Nint)
        s_pass = np.sort([s_pass1,s_pass2])
        Nint = 2 * Nint
        oldT0 = []
        for k in range(0,len(s_pass)):
            Ainput = np.diag(A)
            Einput = np.zeros(Nc)
            Y = fitcalcABCDE(s_pass[k],Ainput,B,C,D,Einput)
            G = np.real(Y)
            if colinterch == 0:
                EE[:,k] = lin.eig(G)
            else:
                [T0,DD] = lin.eig(G)
                T0 = rot(T0)
                [T0,DD] = intercheig(T0,oldT0,DD,Nc,k)
                oldT0 = T0
                EE[:,k] = np.diag(DD)
        
        #===================================================
        #   Indentifying violation, picking minima for s2
        #===================================================
        
        s_pass_ind = np.zeros([1,len(s_pass)])
        for row in range(0,Nc):
            if EE[row,0] < 0:
                s_pass_ind[0] = 1
        for k in range(1,len(s_pass)-1):
            for row in range(0,Nc):
                if EE[row,k] < 0:
                    if(EE[row,k] < EE[row,k-1]) and (EE[row,k]<EE[row,k+1]):
                        s_pass_ind[k] = 1
        
        s = [s,s_pass[np.where(s_pass_ind == 1)]]
        dum = min(EE)
        [g_pass2, ind] = min(dum)
        smin2 = s_pass[ind]
        [g_pass, ind] = min([g_pass,g_pass2])
        dums = [smin,smin2]
        smin = dums[ind]
        
        g_pass = min(g_pass,min(min(EE)))
        
    s_pass = s
    
    return s_pass, g_pass, smin