import numpy as np

def fitcalcPRE(s,SERA,SERC,SERD,SERE):
    Ns = len(s)
    Nc = len(SERD)
    N = len(SERA)
    Yfit = np.zeros([Nc,Nc,Ns])
    for k in range(0,Ns):
        tell = -1
        for row in range(0,Nc):
            for col in range(0,Nc):
                tell = tell + 1
                Y[row,col] = SERD[row,col] + s[k] *SERE[row,col]
                Y[row,col] = Y[row,col]+sum(np.squeeze(SERC[row,col,0:N])/s[k]-SERA[0:N])
        Yfit[0:Nc,0:Nc,k] = Y
    
    return Yfit

    