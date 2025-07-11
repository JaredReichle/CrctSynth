import numpy as np
import numpy.linalg as lin
from math import sqrt
from fitcalcPRE import fitcalcPRE
from pr2ss import pr2ss

class FRPOpts:
    solver = 'QUADPROG'
    test = 'None'
    TOLGD = 1e-6
    TOLE = 1e-12
    weightfactor = 0.001
    weightparam = 1
    outputlevel = 1
    auxflag = 1
    

def FRPY(SER,s,s2,s3,opts):
    auxflag = opts.auxflag
    weightfactor = opts.weightfactor
    weightparam = opts.weightparam
    bigweight = opts.weight
    TOLE = opts.TOLE
    TOL = opts.TOLGD
    
    SERA = SER.poles
    [m,n] = SERA.reshape(-1,1).shape
    if m < n:
        SERA = np.transpose(SERA)
    SERC = SER.R
    SERD = SER.D
    SERE = SER.E
    
    opts.H = []
    opts.oldDflag = -1
    opts.oldEflag = -1
    
    #d = lin.eig(SERD)
    #eigD = d
    [eigD,VD] = lin.eig(SERD)
    if any(eigD < 0):
        Dflag = 1
        #[VD,eigD] = lin.eig(SERD)
        invVD = VD**(-1)
        eigD = lin.diag(eigD)
    else:
        Dflag = 0
        
    #e = lin.eig(SERE)
    #eigE = e
    [eigE, VE] = lin.eig(SERE)
    if any(eigE < 0):
        Eflag = 1
        #[VE,eigE] = lin.eig(SERE)
        invVE = VE**(-1)
        eigE = lin.diag(eigE)
    else:
        Eflag = 0
        
    SERCnew = SERC
    SERDnew = SERD
    SEREnew = SERE
    
    N = len(SERA)
    
    bigB = []
    bigC = []
    s = np.transpose(s)
    Ns = len(s)
    Ns2 = len(s2)
    Nc = len(SERD)
    Nc2 = Nc*Nc
    #I = lin.eye(Nc)
    #Mmat = []
    
    #=============================================
    #   Finding out which poles are complex
    #=============================================
    
    #LS Problem
    cindex = np.zeros([N,1])
    for m in range(0,N):
        if np.imag(SERA[m]) != 0:
            if m == 0:
                cindex[m] = 1
            else:
                if cindex[m-1] == 0 or cindex[m-1] == 2:
                    cindex[m] = 1
                    cindex[m+1] = 2
                else:
                    cindex[m] = 2
    
    #======================================
    #   LOOP FOR LEAST SQUARES PROBLEM
    #======================================
    
    bigV = np.zeros([Nc,N*Nc], dtype = 'complex128')
    bigD = np.zeros([Nc,N], dtype = 'complex128')
    biginvV = np.zeros([Nc,N*Nc], dtype = 'complex128')
    
    if(Dflag+Eflag) == 2:
        bigA = np.zeros([Ns*Nc2,Nc*(N+2)], dtype = 'complex128')
    elif (Dflag+Eflag) == 1:
        bigA = np.zeros([Ns*Nc2,Nc*(N+1)], dtype = 'complex128')
    else:
        bigA = np.zeros([Ns*Nc2,Nc*(N)], dtype = 'complex128')
        
    for m in range(0,N):
        R = (SERC[:,:,m])
        if cindex[m] == 0:  #Real pole
            R = R
        elif cindex[m] == 1:    #complex pole, 1st part
            R = np.real(R)
        else:
            R = np.imag(R)
            
        [V,D] = lin.eig(R)
        bigV[0:Nc,m*Nc:(m+1)*Nc] = V
        biginvV[0:Nc,m*Nc:(m+1)*Nc] = V**(-1)
        bigD[:,m] = np.diag(D)
        
    for k in range(0,Ns):
        sk = s[k]
        #Calculating matrix Mmat
        tell = 0
        offs = 0
        Yfit = fitcalcPRE(sk,SERA,SERC,SERD,SERE)
        
        if weightparam != 0:
            if weightparam == 1:
                weight = np.ones([Nc,Nc])
            elif weightparam == 2:
                weight = 1/abs(Yfit)
            elif weightparam == 3:
                weight = 1/sqrt(abs(Yfit))
            elif weightparam == 4:
                weight = np.ones([Nc,Nc])/lin.norm(abs(Yfit))
            elif weightparam == 5:
                weight = np.ones([Nc,Nc])/sqrt(lin.norm(abs(Yfit)))
            else:
                print('ERROR!')
                weight = bigweight[:,:,k]
        else:
            weight = bigweight[:,:,k]
        
        Mmat = np.zeros([Nc*Nc,(Nc*N)], dtype = 'complex128')
        
        for m in range(0,N):
            #V = bigV[:,m*Nc:(m+1)*Nc]
            V = np.squeeze(bigV[:,m*Nc:(m+1)*Nc])
            invV = V**(-1)
            if cindex[m] == 0:  #real pole
                dum = 1/(sk-SERA[m])
            elif cindex[m] == 1:    #complex pole, 1st part
                dum = (1/(sk-SERA[m])) + 1/(sk-np.transpose(SERA[m]))
            else:
                dum = (1j/(sk-np.transpose(SERA[m]))) - 1j/(sk-SERA[m])
            
            for egenverdi in range(0,Nc):
                tell = -1
                gamm = V[:,egenverdi].reshape(-1,1)*invV[egenverdi,:].reshape(1,-1)
                for row in range(0,Nc):
                    for col in range(0,Nc):
                        faktor = weight[row,col]
                        tell = tell + 1
                        if cindex[m] == 0:  #real pole
                            Mmat[tell,offs+egenverdi] = gamm[row,col]*faktor*dum
                        elif cindex[m] == 1:    #complex pole
                            Mmat[tell,offs+egenverdi] = gamm[row,col]*faktor*dum
                        else:
                            Mmat[tell,offs+egenverdi] = gamm[row,col]*faktor*dum
            
            offs = offs + Nc
        
        if Dflag == 1:
            for egenverdi in range(0,Nc):
                gamm = VD[:,egenverdi].reshape(-1,1)*invVD[egenverdi,:].reshape(1,-1)
                tell = -1
                for row in range(0,Nc):
                    for col in range(0,Nc):
                        tell = tell + 1
                        faktor = weight[row,col]
                        Mmat[tell,offs+egenverdi] = gamm[row,col]*faktor
        
        if Eflag == 1:
            for egenverdi in range(0,Nc):
                gamm = VE[:,egenverdi].reshape(-1,1)*invVE[egenverdi,:].reshape(1,-1)
                tell = -1
                for row in range(0,Nc):
                    for col in range(0,Nc):
                        tell = tell + 1
                        faktor = weight[row,col]
                        Mmat[tell,offs+Nc*Dflag+egenverdi] = gamm[row,col]*sk*faktor
        
        bigA[k*Nc2:(k+1)*Nc2,:] = Mmat
        
    #==============================================================
    #   INTRODUCING SAMPLES OUTSIDE LS REGION: ONE SAMPLE PER POLE
    #==============================================================
    
    if auxflag == 1:
        s4 = np.zeros(len(SERA))
        tell = -1
        for m in range(0,len(SERA)):
            if cindex[m] == 0:
                if((abs(SERA[m])>s[Ns-1]/1j)) or (abs(SERA[m]) < s[0]/1j):
                    tell = tell + 1
                    s4[tell] = 1j*abs(SERA[m])
            elif cindex[m] == 1:    #complex pole
                if((abs(np.imag(SERA[m])))>s[Ns-1]/1j) or (abs(np.imag(SERA[m]))<s[0]/1j):
                    tell = tell + 1
                    s4[tell] = 1j*abs(np.imag(SERA[m]))
    
        Ns4 = len(s4)
        
        bigA2 = np.zeros([Ns4*Nc2,Nc*(N+Dflag+Eflag)], dtype = 'complex128')
        
        for k in range(0,Ns4):
            sk = s4[k]
            #Calc matrix Mmat
            tell = -1
            offs = 0
            Yfit = fitcalcPRE(sk,SERA,SERC,SERD,SERE)
            if weightparam == 1:
                weight = np.ones([Nc,Nc])
            elif weightparam == 2:
                weight = 1/abs(Yfit)
            elif weightparam == 3:
                weight = 1/sqrt(abs(Yfit))
            elif weightparam == 4:
                weight = np.ones([Nc,Nc])/lin.norm(abs(Yfit))
            elif weightparam == 5:
                weight = np.ones([Nc,Nc])/sqrt(lin.norm(abs(Yfit)))
            else:
                print('ERROR')
                weight = 1
            weight = weight*weightfactor
            
            for m in range(0,N):
                V = np.squeeze(bigV[:,m*Nc:(m+1)*Nc])
                invV = V**(-1)
                if cindex[m] == 0:  #Real pole
                    dum = 1/(sk-SERA[m])
                elif cindex[m] == 1:
                    dum = (1/(sk-SERA[m]) + 1/(sk-np.transpose(SERA[m])))
                else:
                    dum = (1j/(sk-np.transpose(SERA[m]))) - 1j/(sk-SERA[m])
                
                for egenverdi in range(0,Nc):
                    tell = -1
                    gamm = V[:,egenverdi].reshape(-1,1)*invV[egenverdi,:].reshape(1,-11)
                    for row in range(0,Nc):
                        for col in range(0,Nc):
                            faktor = weight[row,col]
                            tell = tell + 1
                            if cindex[m] == 0:
                                Mmat[tell,offs+egenverdi] = gamm[row,col]*faktor*dum
                            elif cindex[m] == 1:
                                Mmat[tell,offs+egenverdi] = gamm[row,col]*faktor*dum
                            else:
                                Mmat[tell,offs+egenverdi] = gamm[row,col]*faktor*dum
                
                offs = offs + Nc
                
                if Dflag == 1:
                    for egenverdi in range(0,Nc):
                        gamm = VD[:,egenverdi].reshape(-1,1)*invVD[egenverdi,:].reshape(1,-1)
                        tell = -1
                        for row in range(0,Nc):
                            for col in range(0,Nc):
                                tell = tell + 1
                                faktor = weight[row,col]
                                Mmat[tell,offs+egenverdi] = gamm[row,col]*faktor
                                
                if Eflag == 1:
                    for egenverdi in range(0,Nc):
                        gamm = VE[:,egenverdi].reshape(-1,1)*invVE[egenverdi,:].reshape(1,-11)
                        tell = -1
                        for row in range(0,Nc):
                            for col in range(0,Nc):
                                tell = tell + 1
                                faktor = weight[row,col]
                                Mmat[tell, offs+Nc*Dflag+egenverdi] = gamm[row,col]*sk*faktor
                
                bigA2[k*Nc2:(k+1)*Nc2,:] = Mmat
            bigA = np.concatenate([bigA,bigA2])
    
  
    bigA = np.concatenate([np.real(bigA),np.imag(bigA)])
    Acol = len(bigA[0,:])
    Escale = np.zeros(Acol)  
    for col in range(0,Acol):
        Escale[col] = lin.norm(bigA[:,col],2)
        bigA[:,col] = bigA[:,col]/Escale[col]
    bigAdum = np.transpose(bigA)
    H = np.matmul(bigAdum,bigA) #np.transpose(bigA)*bigA
    opts.H = H
    opts.Escale = Escale
    opts.bigV = bigV
    opts.biginvV = biginvV
    if opts.outputlevel == 1:
        print('Done')
    else:
        bigV = opts.bigV
        biginvV = opts.biginvV
        if Dflag != opts.oldDflag or Eflag != opts.oldEflag:
            opts.H = opts.H[0:Nc*(N+Dflag+Eflag),0:Nc*(N+Dflag+Eflag)]
            opts.Escale = opts.Escale[0:Nc*(N+Dflag+Eflag)]
    
    Mmat2 = np.zeros([Nc2,Nc*(N+Dflag+Eflag)], dtype = 'complex128')
    viol_G = []
    viol_D = []
    viol_E = []
    
    #=============================================
    #   LOOP FOR CONSTRAINT PROBLEM, TYPE #1
    #=============================================
    Y = np.zeros([Nc, Ns2], dtype = 'complex128')
    EE = []
    Q = []
    for k in range(0,Ns2):
        if (len(s2)-1) == 0:    #Hard coded - out of loop if no violating intervals
            break
        sk = s2[k]
        for row in range(0,Nc):
            for col in range(0,Nc):
                Y[row,col] = SERD[row,col] + sk*SERE[row,col]
                Y[row,col] = Y[row,col] + sum(np.squeeze(SERC[row,col,0:N])/sk-SERA[0:N])
        
        #Calculating eigenvalues and eigenvectors
        [V,Z] = lin.eig(np.real(Y))
        Z = lin.diag(Z)
        EE[:,k] = np.real(Z)
        if min(np.real(Z)) < 0:
            #calculating matrix M2mat
            tell = -1
            offs = 0
            
            for m in range(0,N):
                VV = bigV[:,(m-1)*Nc+1:m*Nc]
                invVV = biginvV[:,(m-1)*Nc+1:m*Nc]
                for egenverdi in range(0,Nc):
                    tell = -1
                    gamm = VV[:,egenverdi]*invVV[egenverdi,:]
                    for row in range(0,Nc):
                        for col in range(0,Nc):
                            tell = tell + 1
                            if cindex[m] == 0:
                                Mmat2[tell,offs+egenverdi] = gamm[row,col]/(sk-SERA[m])
                            elif cindex[m] == 1:
                                Mmat2[tell,offs+egenverdi] = gamm[row,col]*(1/(sk-SERA[m]) + 1/(sk-np.transpose(SERA[m])))
                            else:
                                Mmat2[tell,offs+egenverdi] = gamm[row,col]*(1j/(sk-np.transpose(SERA[m])) - 1j/(sk-SERA[m]))
                
                offs = offs + Nc
            
            if Dflag == 1:
                for egenverdi in range(0,Nc):
                    tell = -1
                    gamm = VD[:,egenverdi]*invVD[egenverdi,:]
                    for row in range(0,Nc):
                        for col in range(0,Nc):
                            tell = tell + 1
                            Mmat2[tell,offs+egenverdi] = gamm[row,col]
            
            for n in range(0,Nc):
                tell = -1
                V1 = V[:,n]
                for row in range(0,Nc):
                    for col in range(0,Nc):
                        if row == col:
                            qij = V1[row]**2
                        else:
                            qij = V1[row]*V1[col]
                        tell = tell + 1
                        Q[n,tell] = qij
            
            B = Q*Mmat2
            delz = np.real(Z)
            for n in range(0,Nc):
                if delz[n] > 0:
                    bigB = [bigB,B[n,:]]
                    bigC = [bigC,-TOL+delz[n]]
                    viol_G = [viol_G, delz[n]]
            
    
    #===============================
    #   LOOP FOR CONSTRAINT PROBLEM
    #===============================
    Ns3 = len(s3)
    for k in range(0,Ns3):
        if (Ns3-1) == 0:
            break
        sk = s3[k]
        for row in range(0,Nc):
            for col in range(0,Nc):
                Y[row,col] = SERD[row,col] + sk*SERE[row,col]
                Y[row,col] = Y[row,col] + sum(np.squeeze(SERC[row,col,0:N])/(sk-SERA[0:N]))
        
        #Calculating eigenvalues and eigenvectors
        [V,Z] = lin.eig(np.real(Y))
        Z = lin.diag(Z)
        EE[:,k] = np.real(Z)
        
        tell = -1
        offs = 0
        for m in range(0,N):
            VV = bigV[:,(m-1)*Nc+1:m*Nc]
            invVV = biginvV[:,(m-1)*Nc+1:m*Nc]
            for egenverdi in range(0,Nc):
                tell = -1
                gamm = VV[:,egenverdi]*invVV[egenverdi,:]
                for row in range(0,Nc):
                    for col in range(0,Nc):
                        tell = tell + 1
                        if cindex[m] == 0:
                            Mmat2[tell,offs+egenverdi] = gamm[row,col]/(sk-SERA[m])
                        elif cindex[m] == 1:
                            Mmat2[tell,offs+egenverdi] = gamm[row,col]*((1/(sk-SERA[m]))+ 1/(sk-np.transpose(SERA[m])))
                        else:
                            Mmat2[tell,offs+egenverdi] = gamm[row,col]*(1j/(sk-np.transpose(SERA[m]))+ 1j/(sk-SERA[m]))
            
            offs = offs + Nc
        
        if Dflag == 1:
            for egenverdi in range(0,Nc):
                tell = -1
                gamm = VD[:,egenverdi]*invVD[egenverdi,:]
                for row in range(0,Nc):
                    for col in range(0,Nc):
                        tell = tell + 1
                        Mmat2[tell,offs+egenverdi] = gamm[row,col]
        
        for n in range(0,Nc):
            tell = -1
            V1 = V[:,n]
            for row in range(0,Nc):
                for col in range(0,Nc):
                    if row == col:
                        qij = V1[row]**2
                    else:
                        qij = V1[row]*V1[col]
                    tell = tell + 1
                    Q[n,tell] = qij
        
        B = Q*Mmat2
        delz = np.real(Z)
        for n in range(0,Nc):
            bigB = [bigB,B[n,:]]
            bigC = [bigC,-TOL+delz[n]]
            viol_G = [viol_G, delz[n]]
    
    if Dflag == 1:
        for n in range(0,Nc):
            dum = np.zeros([1,(Nc*(N+Dflag+Eflag))])
            dum[Nc*N+n] = 1
            bigB = [bigB,dum]
            bigC = [bigC,eigD[n]-TOL]
            viol_G = [viol_G, eigD[n]]
            viol_D = [viol_D, eigD[n]]
    
    if Eflag == 1:
        for n in range(0,Nc):
            dum = np.zeros([1,(Nc*(N+Dflag+Eflag))])
            dum[Nc*(N+Dflag)+n] = 1
            bigB = [bigB,dum]
            bigC = [bigC,eigE[n]-TOLE]
            viol_E = [viol_E, eigE[n]]
    if len(bigB) == 0:
        return SER, opts     #No passivity violations
    
    #C = bigC
    
    bigB = [np.real(bigB)]
    for col in range(0,len(opts.H)):
        if len(bigB) > 0:
            bigB[:,col] = bigB[:,col]/opts.Escale[col]
    
    ff = np.zeros([len(opts.H),1])
    
    
    #FIND ROUTINE FOR THE BELOW
    if opts.solver == 'QUADPROG':
        [dx,lambd] = quadprog_solve_qp(opts.H,ff,-bigB,bigC)
    elif opts.solver == 'CPLEX':
        print('ERROR: Python does not support CPLEX solver')
        return
    
    dx = dx/np.transpose(opts.Escale)
    
    for m in range(0,N):
        if cindex[m] == 0:
            D1 = lin.diag(dx[(m-1)*Nc+1:m*Nc])
            SERCnew[:,:,m] = SERCnew[:,:,m] + bigV[:,(m-1)*Nc+1:m*Nc]*D1*biginvV[:,(m-1)*Nc+1:m*Nc]
        elif cindex[m] == 1:
            GAMM1 = bigV[:,(m-1)*Nc+1:m*Nc]
            GAMM2 = bigV[:,(m+1-1)*Nc+1:(m+1)*Nc]
            invGAMM1 = biginvV[:,(m-1)*Nc+1:m*Nc]
            invGAMM2 = biginvV[:,(m+1-1)*Nc+1:(m+1)*Nc]
            
            D1 = lin.diag(dx[(m-1)*Nc+1:m*Nc])
            D2 = lin.diag(dx[(m+1-1)*Nc+1:(m+1)*Nc])
            R1 = np.real(SERC[:,:,m])
            R2 = np.imag(SERC[:,:,m])
            R1new = R1 + GAMM1*D1*invGAMM1
            R2new = R2 + GAMM2*D2*invGAMM2
            SERCnew[:,:,m] = R1new+1j*R2new
            SERCnew[:,:,m+1]
    
    if Dflag == 1:
        DD = lin.diag(dx[N*Nc+1:(N+1)*Nc])
        SERDnew = SERDnew+VD*DD*invVD
    if Eflag == 1:
        EE = lin.diag(dx[N+Dflag*Nc+1:(N+Dflag+Eflag)*Nc])
        SEREnew = SEREnew+VE*EE*invVE
    
    SERDnew = (SERDnew+np.transpose(SERDnew))/2
    SEREnew = (SEREnew+np.transpose(SEREnew))/2
    for m in range(0,N):
        SERCnew[:,:,m] = (SERCnew[:,:,m]+np.transpose(SERCnew[:,:,m]))/2
    
    SER.R = SERCnew
    SER.D = SERDnew
    SER.E = SEREnew
    [SER] = pr2ss(SER)
    
    opts.oldDflag = Dflag
    opts.oldEflag = Eflag
    
    return SER, opts

def quadprog_solve_qp(P, q, G=None, h=None, A=None, b=None):
    qp_G = .5 * (P + P.T)   # make sure P is symmetric
    qp_a = -q
    if A is not None:
        qp_C = -np.vstack([A, G]).T
        qp_b = -np.hstack([b, h])
        meq = A.shape[0]
    else:  # no equality constraint
        qp_C = -G.T
        qp_b = -h
        meq = 0
    return quadprog.solve_qp(qp_G, qp_a, qp_C, qp_b, meq)[0]
            
            