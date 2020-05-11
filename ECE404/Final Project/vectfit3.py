import numpy as np
import numpy.linalg as lin
from math import sqrt, pi
from matplotlib import pyplot as plt

class SERClass:
    def __init__(self, A, B, C, D, E):
        self.A = A
        self.B = B
        self.C = C
        self.D = D
        self.E = E
        self.R = None
        self.poles = None

class DefaultVect3Opts():
    relax=1      #Use vector fitting with relaxed non-triviality constraint
    stable=1     #Enforce stable poles
    asymp=1      #Include only D in fitting (not E)   
    skip_pole=0  #Do NOT skip pole identification
    skip_res=0   #Do NOT skip identification of residues (C,D,E) 
    cmplx_ss=1   #Create complex state space model
    spy1=0       #No plotting for first stage of vector fitting
    spy2=1       #Create magnitude plot for fitting of f(s) 
    logx=1       #Use logarithmic abscissa axis
    logy=1       #Use logarithmic ordinate axis 
    errplot=1    #Include deviation in magnitude plot
    phaseplot=0  #exclude plot of phase angle (in addition to magnitiude)
    legend=1

def vectfit3(f,s,poles,weight,opts):
    
    
    #Tolerances used by relaxed version of vector fitting
    TOLlow = 1e-18
    TOLhigh = 1e18
    [a,b] = poles.shape
    if a > b:
        poles = np.transpose(poles)
    if s[0,0] == 0 and a == 1:
        if poles[0] == 0 and poles[1] != 0:
            poles[0] = -1
        elif poles[1] == 0 and poles[0] != 0:
            poles[1] = -1
        elif poles[0] == 0 and poles[1] == 0:
            poles[0] = -1+1j*10
            poles[1] = -1-1j*10
    
    if (opts.relax != 0) and (opts.relax != 1):
        print('ERROR in vectfit3: ==> Illegal value for opts.relax')
        return
    if (opts.asymp != 1) and (opts.asymp != 2) and (opts.asymp != 3):
        print('ERROR in vectfit3: ==> Illegal value for opts.asymp')
        return
    if (opts.stable != 0) and (opts.stable != 1):
        print('ERROR in vectfit3: ==> Illegal value for opts.stable')
        return
    if (opts.skip_pole != 0) and (opts.skip_pole != 1):
        print('ERROR in vectfit3: ==> Illegal value for opts.skip_pole')
        return
    if (opts.skip_res != 0) and (opts.skip_res != 1):
        print('ERROR in vectfit3: ==> Illegal value for opts.skip_res')
        return
    if (opts.cmplx_ss != 0) and (opts.cmplx_ss != 1):
        print('ERROR in vectfit3: ==> Illegal value for opts.cmplx_ss')
        return

    rmserr = []
    [a,b] = s.shape
    if a < b:
        s = np.transpose(s)
        
    #Sanity checks on input arrays
    if len(s) != len(f[0]):
        print('ERROR in vectfit3: ==> Second dimension of f does not match length of s.')
        return
    if len(s) != len(weight[0,:]):
        print('ERROR in vectfit3: ==> Second dimension of weight does not match length of s')
        return
    if len(weight[:,0]) != 1:
        if len(weight[:,0]) != len(f):
            print('ERROR in vectfit3: ==> First dimension of weight is neither 1 nor matches first dimension of f')
            return
    
    
    LAMBD = np.diag(poles[0]) #NxN diagonal
    Ns = len(s)
    N = len(LAMBD)
    Nc = len(f) #HEY LOOK <--- ITS THIS GUY AGAIN
    B = np.ones([N,1])
    SERA = poles
    SERC = np.zeros([Nc,N])
    SERD = np.zeros(Nc)
    SERE = np.zeros(Nc)
    roetter = poles
    fit = np.zeros([Nc,Ns], dtype = 'complex128')
    
    weight = np.transpose(weight)
    if len(weight[0,:]) == 1:
        common_weight = 1
    elif len(weight[0,:]) == Nc:
        common_weight = 0
    else:
        print('ERROR in vectfit3: Invalid size of array weight')
        return
    
    if opts.asymp == 1:
        offs = 0
    elif opts.asymp == 2:
        offs = 1
    else:
        offs = 2
    
    
    #====================================
    #   POLE IDENTIFICATION
    #====================================
    
    if opts.skip_pole != 1:
        Escale = np.zeros([1,Nc+1])
        
        #==================================================
        #   Finding out which starting poles are complex
        #==================================================
        cindex = np.zeros([N,1])
        for m in range(0,N):
            if np.imag(LAMBD[m,m]) != 0:
                if m == 0:
                    cindex[m] = 1 #First part complex
                else:
                    if cindex[m-1] == 0 or cindex[m-1] == 2: #One before real, or one before second complex
                       cindex[m] = 1 #First part complex
                       if m != N-1:
                           cindex[m+1] = 2 #Second part complex
                    else:   #One before is a first part
                        cindex[m] = 2 #Second part complex
        
        #===========================================
        #   Building system - matrix
        #===========================================
        
        Dk = np.zeros([Ns,N], dtype = 'complex128')
        for m in range(0,N):
            if cindex[m] == 0:  #Real pole
                #tmpreq = 1/(s-LAMBD[m,m])
                #tmpreq = tmpreq.reshape(-1,1)
                Dk[:,m] = (1/(s-LAMBD[m,m]))[:,0]
            elif cindex[m] == 1: #Complex pole, first part
                tmp1 = 1/(s-LAMBD[m,m])
                tmp2 = 1/(s-np.conj(LAMBD[m,m]))
                tmp3 = tmp2 + tmp1
                Dk[:,m] = tmp3[:,0]
                #Dk[:,m] = 1/(s-LAMBD[m,m]) + 1/(s-LAMBD[m,m])
                tmp1 = 1j/(s-LAMBD[m,m])
                tmp2 = 1j/(s-np.conj(LAMBD[m,m]))
                tmp3 = tmp1 - tmp2
                Dk[:,m+1] = tmp3[:,0]
                #Dk[:,m+1] = 1j/(s-LAMBD[m,m]) - 1j/(s-np.conj(LAMBD[m,m]))
        if opts.asymp == 1 or opts.asymp == 2:
            Dk = np.hstack((Dk,np.ones([Ns,1])))
        elif opts.asymp == 3:
            Dk = np.hstack((Dk,np.ones([Ns,1])))
            Dk = np.hstack((Dk,s[:,0]))
        
        #Scaling for last row of LS-problem (pole id)
        scale = 0
        for m in range(0,Nc):
            if len(weight[0,:]) == 1:
                scale = scale+(lin.norm(weight[:,0]*np.transpose(f[m,:])))**2
            else:
                scale = scale+(lin.norm(weight[:,m]*np.transpose(f[m,:])))**2
        scale = sqrt(scale)/Ns
        
        if opts.relax == 1:
            AA = np.zeros([Nc*(N+1),N+1], dtype = 'complex128')
            bb = np.zeros([Nc*(N+1),1])
            Escale = np.zeros([1, len(AA[0,:])])
            for n in range(0,Nc):
                A = np.zeros([Ns,(N+offs)+N+1], dtype = 'complex128')
                if common_weight == 1:
                    weig = weight
                else:
                    weig = weight[:,n]
                
                for m in range(0,N+offs):#Left block
                    tmp1 = weig
                    tmp2 = Dk[:,m].reshape(-1,1)
                    tmp3 = tmp1*tmp2
                    A[:,m] = tmp3[:,0]
                    #A[:,m] = weig[:,0]*Dk[:,m] #302x21 = 302x1 times 302x11
                inda = N+offs
                for m in range(0,N+1): #Right block
                    tmp1 = -weig
                    tmp2 = Dk[:,m].reshape(-1,1)
                    tmp3 = f[n,:].reshape(-1,1)
                    tmp4 = tmp1*tmp2*tmp3
                    A[:,inda+m] = tmp4[:,0]
                    #A[:, inda+m] = weig*Dk[0:Ns,m]*np.transpose(f[n,0:Ns])
                
                A = np.vstack((np.real(A),np.imag(A)))
                
                #Integral criterion for sigma
                offset = (N+offs)
                if (n + 1) == Nc:
                    tmprow = np.zeros([1,A.shape[1]])
                    for mm in range(0,N+1):
                        tmp1 = sum(Dk[:,mm])
                        tmp2 = np.real(scale*tmp1)
                        tmprow[0,offset+mm] = tmp2
                    A = np.vstack((A,tmprow))
                    #A[2*Ns+1,offset+mm] = np.real(scale*sum(Dk[:,mm]))
                [Q,R] = np.linalg.qr(A)
                ind1 = N+offs #10 + 0 + 1 = 10
                ind2 = N+offs+N+1 #10 + 10 + 1 = 20
                R22 = R[ind1:ind2,ind1:ind2]
                AA[(n)*(N+1):(n+1)*(N+1),:] = R22
                if (n + 1) == Nc:
                    tmpq1 = Q[-1,N+offs:]
                    tmpq2 = np.transpose(tmpq1)
                    tmpq3 = Ns*scale*tmpq2
                    tmpindb1 = n*(N+1)
                    tmpindb2 = (n+1)*(N+1)
                    bb[tmpindb1:tmpindb2,0] = tmpq3
                    #bb[n*(N+1):(n+1)*(N+1),0] = Q[-1,N+offs:]
            
            for col in range(0,len(AA[0,:])):
                Escale[0,col] = 1/lin.norm(AA[:,col])
                AA[:,col] = Escale[0,col]*AA[:,col]
            x = lin.solve(AA,bb)
            x = x*Escale.reshape(-1,1)
        
        if opts.relax == 0 or abs(x[-1]) < TOLlow or abs(x[-1]) > TOLhigh:
            AA = np.zeros([Nc*(N),N], dtype = 'complex128')
            bb = np.zeros([Nc*(N),1])
            if opts.relax == 0:
                Dnew = 1
            else:
                if x[-1] == 0:
                    Dnew = 1
                elif abs(x[-1]) < TOLlow:
                    Dnew = np.sign(x[-1])*TOLlow
                elif abs(x[-1]) > TOLhigh:
                    Dnew = np.sign(x[-1])*TOLhigh
            
            for n in range(0,Nc):
                A = np.zeros([Ns,(N+offs)+N], dtype = 'complex128')
                Escale = np.zeros([1,N], dtype = 'complex128')
                if common_weight == 1:
                    weig = weight
                else:
                    weig = weight[:,n]
                    
                for m in range(0,N+offs):
                    tmp1 = weig
                    tmp2 = Dk[:,m].reshape(-1,1)
                    tmp3 = tmp1*tmp2
                    A[:,m] = tmp3[:,0]
                inda = N+offs
                for m in range(0,N):
                    tmp1 = weig
                    tmp2 = Dk[:,m].reshape(-1,1)
                    tmp3 = f[n,:].reshape(-1,1)
                    tmp4 = tmp1*tmp2*tmp3
                    A[:,inda+m] = tmp4[:,0]
                    #A[0:Ns, inda+m] = -weig*Dk[0:Ns,m]*np.transpose(f[n,0:Ns])
                b = Dnew*weig*f[n,:].reshape(-1,1)
                #b = Dnew*weig*np.transpose(f[n,0:Ns])
                A = np.vstack((np.real(A),np.imag(A)))
                b = np.vstack((np.real(b),np.imag(b)))
                offset = (N+offs)
                [Q,R] = np.linalg.qr(A)
                ind1 = N+offs
                ind2 = N+offs+N
                R22 = R[ind1:ind2,ind1:ind2]
                AA[n*N:(n+1)*N,:] = R22
                tmpnew1 = Q[:,ind1:ind2]
                tmpnew2 = np.transpose(tmpnew1)
                tmpnew3 = np.matmul(tmpnew2,b)
                bb[n*N:(n+1)*N,0] = tmpnew3[0]
                #bb[n*N:(n+1)*N,0] = np.transpose(Q[:,ind1:ind2])*b #10x604 times 604x1
            for col in range(0,len(AA[0,:])):
                Escale[0,col] = 1/lin.norm(AA[:,col])
                AA[:,col] = Escale[0,col]*AA[:,col]
            
            x = lin.solve(AA,bb)
            x = x*Escale.reshape(-1,1)
            x = np.vstack((x,Dnew))
            
            
        C = x[0:-1]
        D = x[-1]
        
        RES3 = []
        
        #We now change back to make C complex:
        for m in range(0,N):
            if cindex[m] == 1:
                r1 = np.real(C[m])
                r2 = np.real(C[m+1])
                C[m] = r1+1j*r2
                C[m+1] = r1-1j*r2
        
        if opts.spy1 == 1:
            Dk = np.zeros([Ns,N])
            for m in range(0,N):
                Dk[:,m] = 1/(s-LAMBD[m,m])
            RES3[:,0] = D + Dk*C
            
            freq = s/(2*pi*1j)
            
            if opts.logx == 1:
                if opts.logy == 1:
                    plt.figure(3)
                    plt.loglog(freq,abs(np.conj(np.transpose(RES3))), 'b')
                    plt.xlim([freq[0], freq[Ns]])
                else: #logy = 0
                    plt.figure(3)
                    plt.semilogx(freq, abs(np.conj(np.tranpose(RES3))), 'b')
                    plt.xlim([freq[0], freq[Ns]])
            else: #logx = 0
                if opts.logy == 1:
                    plt.figure(3)
                    plt.semilogy(freq,abs(np.conj(np.transpose(RES3))),'b')
                    plt.xlim([freq[0],freq[Ns]])
                else:
                    plt.figure(3)
                    plt.plot(s/(2*pi*1j),abs(np.conj(np.transpose(RES3))),'b')
                    plt.xlim([freq[0],freq[Ns]])
            
            plt.xlabel('Frequency [Hz]')
            plt.ylabel('Magnitude')
            if opts.legend == 1:
                plt.legend('sigma')
            
        #==========================================
        #   We now calculate the zeros for sigma
        #==========================================
        
        m = -1
        
        
        for n in range(0,N):
            m = m + 1
            if m < N:
                if(abs(LAMBD[m,m]) > abs(np.real(LAMBD[m,m]))): #Complex number?
                    LAMBD[m+1,m] = -np.imag(LAMBD[m,m])
                    LAMBD[m,m+1] = np.imag(LAMBD[m,m])
                    LAMBD[m,m] = np.real(LAMBD[m,m])
                    LAMBD[m+1,m+1] = LAMBD[m,m]
                    B[m] = 2
                    B[m+1] = 0
                    koko = C[m]
                    C[m+1] = koko.imag
                    C[m] = koko.real
                    m = m + 1
        
        tmp1 = np.transpose(C) #CHECKED
        tmp2 = B*tmp1      
        tmp3 = tmp2/D[0]    
        tmp4 = LAMBD - tmp3
        ZER = LAMBD - (B*np.transpose(C))/D[0]
        reig, reig2 = lin.eig(ZER)
        roetter = np.transpose(reig)
        unstables = np.real(roetter) > 0
        if opts.stable == 1:
            roetter[unstables] = roetter[unstables] - 2*np.real(roetter[unstables])
        roetter = np.sort(roetter)
        N = len(roetter)
        
        #====================
        for n in range(0,N):
            for m in range(n+1,N):
                if roetter[m].imag == 0 and roetter[n].imag != 0:
                    trans = roetter[n]
                    roetter[n] = roetter[m]
                    roetter[m] = trans
        
        N1 = 0
        for m in range(0,N):
            if roetter[m].imag == 0:
                N1 = m
        if N1 < N:
            roetter[N1:N] = np.sort(roetter[N1:N])
        
        roetter = roetter - 2*1j*np.imag(roetter)
        SERA = np.transpose(roetter)
        
    #=========================================
    #   RESIDUE IDENTIFICATION
    #=========================================
    if opts.skip_res != 1:
        
        #We now calculate SER for f, using the modified zeros of sigma as new poles:
        
        tmproe = roetter.reshape(-1,1)
        LAMBD = tmproe
        
        cindex = np.zeros([N,1])
        for m in range(0,N):
            if np.imag(LAMBD[m]) != 0:
                if m == 0:
                    cindex[m] = 1
                else:
                    if cindex[m-1] == 0 or cindex[m-1] == 2: #One before real, or one before second complex
                       cindex[m] = 1 #First part complex
                       if m != N-1:
                           cindex[m+1] = 2 #Second part complex
                    else:   #One before is a first part
                        cindex[m] = 2 #Second part complex
                        
                        
        
        #====================================================================
        # We now calculate the SER for f (new fitting), using the above calc
        # zeros as known poles :
        #====================================================================
        
        if opts.asymp == 1:
            A = np.zeros([2*Ns,N], dtype = 'complex128') #[2*302,10]
            BB = np.zeros([2*Ns,Nc], dtype = 'complex128') #[2*302,1]
        elif opts.asymp == 2:
            A = np.zeros([2*Ns,N+1], dtype = 'complex128')
            BB = np.zeros([2*Ns,Nc], dtype = 'complex128')
        else:
            A = np.zeros([2*Ns,N+2], dtype = 'complex128')
            BB = np.zeros([2*Ns,Nc], dtype = 'complex128')
        
        Dk = np.zeros([Ns,N], dtype = 'complex128')
        for m in range(0,N):
            if cindex[m] == 0:  #Real pole                
                Dk[:,m] = 1/(s-LAMBD[m])
            elif cindex[m] == 1:    #Complex pole, first part
                tmp1 = weight/(s-LAMBD[m])
                tmp2 = weight/(s-np.conj(LAMBD[m]))
                tmp3 = tmp2 + tmp1
                Dk[:,m] = tmp3[:,0]
                #Dk[:,m] = weight/(s-LAMBD[m]) + weight/(s-np.conj(LAMBD[m]))
                tmp1 = 1j*weight/(s-LAMBD[m])
                tmp2 = 1j*weight/(s-np.conj(LAMBD[m]))
                tmp3 = tmp1 - tmp2
                Dk[:,m] = tmp3[:,0]
                #Dk[:,m+1] = 1j*weight/(s-LAMBD[m]) - 1j*weight/(s-np.conj(LAMBD[m]))
        
        if common_weight == 1:
            
            Dk = np.zeros([Ns,N], dtype = 'complex128')
            for m in range(0,N):
                if cindex[m] == 0:  #Real pole
                    Dk[:,m] = weight/(s-LAMBD[m])
                elif cindex[m] == 1:    #Complex pole, first part
                    tmp1 = weight/(s-LAMBD[m])
                    tmp2 = weight/(s-np.conj(LAMBD[m]))
                    tmp3 = tmp2 + tmp1
                    Dk[:,m] = tmp3[:,0] #ISSUE - Complex casting
                    #Dk[:,m] = weight/(s-LAMBD[m]) + weight/(s-np.conj(LAMBD[m]))
                    tmp1 = 1j*weight/(s-LAMBD[m])
                    tmp2 = 1j*weight/(s-np.conj(LAMBD[m]))
                    tmp3 = tmp1 - tmp2
                    Dk[:,m+1] = tmp3[:,0] #ISSUE - Complex casting
                    #Dk[:,m+1] = 1j*weight/(s-LAMBD[m]) - 1j*weight/(s-np.conj(LAMBD[m]))
            
            #LINE 600
            
            if opts.asymp == 1:
                A[0:Ns,0:N] = Dk
            elif opts.asymp == 2:
                A[0:Ns,0:N] = Dk
                A[0:Ns,N+1] = weight
            else:
                A[0:Ns,0:N] = Dk
                A[0:Ns,N+1] = weight
                A[0:Ns,N+2] = weight*s
            
            for m in range(0,Nc):
                tmp1 = f[m,:]
                tmp2 = tmp1.reshape(-1,1)
                tmp3 = weight
                tmp4 = tmp3*tmp2
                BB[0:Ns,m] = tmp4[:,0]
                #BB[0:Ns,m] = weight*(f[m,:].reshape(-1,1))
            A[Ns:2*Ns,:] = np.imag(A[0:Ns,:])
            A[0:Ns,:] = np.real(A[0:Ns,:])
            BB[Ns:2*Ns,:] = np.imag(BB[0:Ns,:])
            BB[0:Ns,:] = np.real(BB[0:Ns,:])
            
            if opts.asymp == 2:
                A[0:Ns,N+1] = A[0:Ns,N+1]
            elif opts.asymp == 3:
                A[0:Ns,N+1] = A[0:Ns,N+1]
                A[Ns:2*Ns] = A[Ns:2*Ns,N+2]
                
            #clear Escale
            Escale = np.zeros([len(A[0,:]),1])
            for col in range(0,len(A[0,:])):
                Escale[col] = lin.norm(A[:,col],2)
                A[:,col] = A[:,col]/Escale[col]   
            [X, residuals, rank, sva] = lin.lstsq(A,BB, rcond=None)
            for n in range(0,Nc):
                X[:,n] = X[:,n]/np.transpose(Escale)
            
            #clear A
            X = np.transpose(X)
            C = X[:,0:N]
            if opts.asymp == 2:
                SERD = X[:,N+1]
            elif opts.asymp == 3:
                SERE = X[:,N+2]
                SERD = X[:,N+1]
            
        else:   #if common_weight == 1
            
            SERD = np.zeros(Nc)
            SERE = np.zeros(Nc)
            C = np.zeros([Nc,N])
            for n in range(0,Nc):
                if opts.asymp == 1:
                    A[0:Ns,0:N] = Dk
                elif opts.asymp == 2:
                    A[0:Ns,0:N] = Dk
                    A[0:Ns,N+1] = 1
                else:
                    A[0:Ns,0:N] = Dk
                    A[0:Ns,N+1] = 1
                    A[0:Ns,N+2] = s
                for m in range(0,len(A[0,:])):
                    A[0:Ns,m] = weight[:,n]*A[0:Ns,m]
                
                BB = weight[:,n]*np.transpose(f[n,:])
                A[Ns+1:2*Ns,:] = np.imag(A[0:Ns,:])
                A[0:Ns,:] = np.real(A[0:Ns,:])
                BB[Ns+1:2*Ns] = np.imag(BB[0:Ns])
                BB[0:Ns] = np.real(BB[0:Ns])
                
                if opts.asymp == 2:
                    A[0:Ns,N+1] = A[0:Ns,N+1]
                elif opts.asymp == 3:
                    A[0:Ns,N+1] = A[0:Ns,N+1]
                    A[Ns+1:2*Ns,N+2] = A[Ns+1:2*Ns,N+2]
                
                #clear Escale
                Escale = np.zeros([1,len(A[0,:])])
                for col in range(0,len(A[0,:])):
                    Escale[col] = lin.norm(A[:,col],2)
                    A[:,col] = A[:,col]/Escale[col]
                x = A/BB
                x = x/np.transpose(Escale)
                
                #clear A
                C[n,0:N] = np.transpose(x[0:N])
                
                if opts.asymp == 2:
                    SERD[n] = x[N+1]
                elif opts.asymp == 3:
                    SERE[n] = x[N+2]
                    SERD[n] = x[N+1]
                
            #end for n=0:Nc
        #end if common_weight == 1
        
        #===================================================
        
        #We now change back to make C complex
        for m in range(0,N):
            if cindex[m] == 1:
                for n in range(0,Nc):
                    r1 = C[n,m]
                    r2 = C[n,m+1]
                    C[n,m] = r1+1j*r2
                    C[n,m+1] = r1-1j*r2
        
        
        B = np.ones([N,1])
        
        #===================================================
        
        SERA = LAMBD
        SERB = B
        SERC = C
        
        Dk = np.zeros([Ns,N], dtype = 'complex128')
        for m in range(0,N):
            tmp1 = SERA[m]
            tmp2 = s - tmp1
            tmp3 = 1/tmp2
            Dk[:,m] = tmp3[:,0]
            #Dk[:,m] = 1/(s-SERA[m]) #[302,10] = 1/ ([10,1] - [10,1])
        for n in range(0,Nc):
            tmp1 = SERC[n,:].reshape(1,-1) #Row vector
            tmp2 = np.transpose(tmp1) #Column vector
            tmp3 = np.matmul(Dk,tmp2)
            tmp4 = np.transpose(tmp3)
            fit[n,:] = tmp4
            #fit[n,:] = np.transpose(Dk*np.transpose(SERC[n,:]))
            if opts.asymp == 2:
                fit[n,:] = fit[n,:] + SERD[n]
            elif opts.asymp == 3:
                fit[n,:] = fit[n,:] + SERD[n]+np.tranpose(s)*SERE[n]
        
        #fit = np.transpose(fit)
        diff = fit - f
        rmserr = sqrt(sum(sum(abs(diff**2))))/sqrt(Nc*Ns)
        
        
        #PLOTTING STUFF HERE
        
        fit = np.transpose(fit)
        
    #end if skip_res != 1
    
    A = SERA
    poles = A
    
    if opts.skip_res != 1:
        B = SERB
        C = SERC
        D = SERD
        E = SERE
    else:
        B = np.ones([N,1])
        C = np.zeros([Nc,N])
        D = np.zeros([Nc,Nc])
        E = np.zeros([Nc,Nc])
        rmserr = 0
    
    #=========================================
    #   Convert into real state-space model
    #=========================================
    
    if opts.cmplx_ss != 1:
        A = lin.diag(lin.sparse(A))
        
        cindex = np.zeros([N,1])
        for m in range(0,N):
            if np.imag(A[m,m]) != 0:
                if m == 0:
                    cindex[m] = 1
                else:
                    if cindex[m-1] == 0 or cindex[m-1] == 2:
                        cindex[m] = 1
                        cindex[m+1] = 2
                    else:
                        cindex[m] = 2
        
        n = -1
        for m in range(0,N):
            n = n + 1
            if cindex[m] == 1:
                a = A[n,n]
                a1 = np.real(a)
                a2 = np.imag(a)
                c = C[:,n]
                c1 = np.real(c)
                c2 = np.imag(c)
                b = B[n,:]
                b1 = 2*np.real(b)
                b2 = 2*np.imag(b)
                Ablock = [[a1,a2],[-a2,a1]]
                
                A[n:n+1,n:n+1] = Ablock
                C[:,n] = c1
                C[:,n+1] = c2
                B[n,:] = b1
                B[n+1,:] = b2
        else:
            A = lin.sparse(lin.diag(A))
        
    #end if cmplx_ss != 1
            
    SER = SERClass(A,B,C,D,E)
    
    poles = poles.reshape(-1,1)
        
    return SER, poles, rmserr, fit