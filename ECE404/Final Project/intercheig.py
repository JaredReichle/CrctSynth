import numpy as np

def intercheig(V,oldV,D,Nc,fstep):
    dot = np.zeros(Nc)
    ind = np.zeros(Nc)
    taken = np.zeros(Nc)
    hjelp = np.zeros(Nc)
    if fstep > 0:
        UGH = abs(np.real(np.transpose(oldV)*V))
        for ii in range(0,Nc):
            ilargest = 0
            rlargest = 0
            for j in range(0,Nc):
                dotprod = UGH[ii,j]
                if dotprod > rlargest:
                    rlargest = abs(np.real(dotprod))
                    ilargest = j
            dot[ii] = rlargest
            ind[ii] = ii
            taken[ii] = 0
        
        for ii in range(0,Nc):
            for j in range(0,Nc-1):
                if dot[j] < dot[j+1]:
                    hjelp[0] = dot[j+1]
                    ihjelp = ind[j+1]
                    dot[j+1] = dot[j]
                    ind[j+1] = ind[j]
                    dot[j] = hjelp[0]
                    ind[j] = ihjelp
        
        for l in range(0,Nc):
            ii = int(ind[l])
            ilargest = 0
            rlargest = 0
            
            for j in range(0,Nc):
                if taken[j] == 0:
                    dotprod = UGH[ii,j]
                    if dotprod > rlargest:
                        rlargest = abs(np.real(dotprod))
                        ilargest = j
            
            taken[ii] = 1
            
            hjelp = V[:,ii]
            V[:,ii] = V[:,ilargest]
            V[:,ilargest] = hjelp
            
            hjelp = D[ii]
            D[ii] = D[ilargest]
            D[ilargest] = hjelp
            
            dum = UGH[:,ii]
            UGH[:,ii] = UGH[:,ilargest]
            UGH[:,ilargest] = dum
            
        for ii in range(0,Nc):
            dotprod = np.transpose(oldV[:,ii])*V[:,j]
            if np.sign(np.real(dotprod[ii])) < 0:
                V[:,ii] = -V[:,ii]
        
    return V, D