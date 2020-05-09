import numpy as np
import numpy.linalg as lin
from math import sqrt

class defa():
    relax=1      #Use vector fitting with relaxed non-triviality constraint
    stable=1     #Enforce stable poles
    asymp=2      #Include only D in fitting (not E)   
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
    [a,b] = np.size(poles)
    if s[0] == 0 and a == 1:
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
    [a,b] = np.size(s)
    if a < b:
        s = np.transpose(s)
        
    #Sanity checks on input arrays
    if len(s) != len(f[0,:]):
        print('ERROR in vectfit3: ==> Second dimension of f does not match length of s.')
        return
    if len(s) != len(weight[0,:]):
        print('ERROR in vectfit3: ==> Second dimension of weight does not match length of s')
        return
    if len(weight[:,0]) != 1:
        if len(weight[:,0]) != len(f[:,0]):
            print('ERROR in vectfit3: ==> First dimension of weight is neither 1 nor matches first dimension of f')
            return
    
    LAMBD = lin.diag(poles)
    Ns = len(s)
    N = len(LAMBD)
    Nc = len(f[:,0])
    B = np.ones(N)
    SERA = poles
    SERC = np.zeros(Nc,N)
    SERD = np.zeros(Nc)
    SERE = np.zeros(Nc)
    roetter = poles
    fit = np.zeros(Nc,Ns)
    
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
        cindex = np.zeors([1,N])
        for m in range(0,N):
            if np.imag(LAMBD[m,m]) != 0:
                if m == 0:
                    cindex[m] = 0
                else:
                    if cindex[m-1] == 0 or cindex[m-1] == 2:
                       cindex[m] = 1
                       cindex[m+1] = 2
                    else:
                        cindex[m] = 2
        
        #===========================================
        #   Building system - matrix
        #===========================================
        
        Dk = np.zeros(Ns,N)
        for m in range(0,N):
            if cindex[m] == 0:  #Real pole
                Dk[:,m] = 1/(s-LAMBD[m,m])
            elif cindex[m] == 1: #Complex pole, first part
                Dk[:,m] = 1/(s-LAMBD[m,m]) + 1/(s-np.transpose(LAMBD[m,m]))
                Dk[:,m+1] = 1j/(s-LAMBD[m,m]) - 1j/(s-np.transpose(LAMBD[m,m]))
        if opts.asymp == 1 or opts.asymp == 2:
            Dk[:,N+1] = 1
        elif opts.asymp == 3:
            Dk[:,N+1] = 1
            Dk[:,N+1] = s
        
        #Scaling for last row of LS-problem (pole id)
        scale = 0
        for m in range(0,Nc):
            if len(weight[0,:]) == 1:
                scale = scale+(lin.norm(weight[:,0]*np.transpose(f[m,:])))**2
            else:
                scale = scale+(lin.norm(weight[:,m]*np.transpose(f[m,:])))**2
        scale = sqrt(scale)/Ns
        
        if opts.relax == 1:
            AA = np.zeros(Nc*(N+1),N+1)
            bb = np.zeros(Nc*(N+1),1)
            Escale = np.zeros(1, len(AA[0,:]))
            for n in range(0,Nc):
                A = np.zeros(Ns,(N+offs)+N+1)
                if common_weight == 1:
                    weig = weight
                else:
                    weig = weight[:,n]
                
                for m in range(0,N+offs): #Left block
                    A[0:Ns,m] = weig*Dk[0:Ns,m]
                inda = N+offs
                for m in range(0,N+1):
                    A[0:Ns, inda+m] = weig*Dk[0:Ns,m]*np.transpose(f[n,0:Ns])
                
                A = [np.real(A), np.imag(A)]
                
                #Integral criterion for sigma
                offset = (N+offs)
                if n == Nc:
                    for mm in range(0,N+1):
                        A[2*Ns+1,offset+mm] = np.real(scale*sum(Dk[:,mm]))
                [Q,R] = np.linalg.qr(A,0)
                ind1 = N+offs+1
                ind2 = N+offs+N+1
                R22 = R[ind1:ind2,ind1:ind2]    #FIND PYTHON ROUTINE
                AA[(n-1)*(N+1)+1:n*(N+1),:] = R22
                if n == Nc:
                    bb[(n-1)*(N+1)+1:n*(N+1),0] = Q[-1,N+offs+1:-1] #Investigate
            
            for col in range(0,len(AA[0,:])):
                Escale[col] = 1/lin.norm(AA[:,col])
                AA[:,col] = Escale[col]*AA[:,col]
            x = AA/bb
            x = x*Escale
        
        if opts.relax == 0 or abs(x[-1]) < TOLlow or abs(x[-1])>TOLhigh:
            AA = np.zeros(Nc*(N),N)
            bb = np.zeros(Nc*(N),1)
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
                A = np.zeros(Ns,(N+offs)+N)
                Escale = np.zeros([1,N])
                if common_weight == 1:
                    weig = weight
                else:
                    weig = weight[:,n]
                    
                for m in range(0,N+offs):
                    A[0:Ns,m] = weig*Dk[0:Ns,m]
                inda = N+offs
                for m in range(0,N):
                    A[0:Ns, inda+m] = -weig*Dk[0:Ns,m]*np.transpose(f[n,0:Ns])
                b = Dnew*weig*np.transpose(f[n,0:Ns])
                A = [[np.real(A)],[np.imag(A)]]
                b = [[np.real(b)],[np.imag(b)]]
                offset = (N+offs)
                [Q,R] = np.linalg.qr(A,0)
                ind1 = N+offs+1
                ind2 = N+offs+N
                R22 = R[ind1:ind2,ind1:ind2]
                AA[(n-1)*N+1:n*N,:] = R22
                bb[(n-1)*N+1:n*N,0] = np.transpose(Q[:,ind1:ind2])*b
            for col in range(0,len(AA[0,:])):
                Escale[col] = 1/lin.norm(AA[:,col])
                AA[:,col] = Escale[col]*AA[:,col]
            
            x = AA/bb
            x = x*np.transpose(Escale)
            x = [[x],[Dnew]]
            
        C = x[0:-2]
        D = x[-1]
        RES3 = []
        
        #We now change back to make C complex:
        for m in range(0,N):
            if cindex[m] == 1:
                for n in range(0,0):
                    r1 = C[m]
                    r2 = C[m+1]
                    C[m] = r1+1j*r2
                    C[m+1] = r1-1j*r2
        
        if opts.spy1 == 1:
            Dk = np.zeros(Ns, N)
            for m in range(0,N):
                Dk[:,m] = 1/(s-LAMBD[m,m])
            RES3[:,0] = D + Dk*C
            #freq = s/(2*pi*1j)
            
            ##PLOTTING
            
        #==========================================
        #   We now calculate the zeros for sigma
        #==========================================
        
        m = 0
        
        for n in range(0,N):
            m = m + 1
            if m < N:
                if(abs(LAMBD[m,m])>abs(np.real(LAMBD[m,m]))): #Complex number?
                    LAMBD[m+1,m] = -np.imag(LAMBD[m,m])
                    LAMBD[m,m+1] = np.imag(LAMBD[m,m])
                    LAMBD[m,m] = np.real(LAMBD[m,m])
                    LAMBD[m+1,m+1] = LAMBD[m,m]
                    B[m,0] = 2
                    B[m+1,0] = 0
                    koko = C[m]
                    C[m] = np.real(koko)
                    C[m] = np.imag(koko)
                    m = m + 1
        
        ZER = LAMBD-B*np.transpose(C)/D
        roetter = np.transpose(lin.eig(ZER))
        unstables = np.real(roetter) > 0
        if opts.stable == 1:
            roetter[unstables] = roetter[unstables] - 2*np.real(roetter[unstables])
        roetter = np.sort(roetter)
        N = len(roetter)
        
        #====================
        for n in range(0,N):
            for m in range(n,N):
                if np.imag(roetter[m]) == 0 and np.imag(roetter[n]) != 0:
                    trans = roetter[n]
                    roetter[n] = roetter[m]
                    roetter[m] = trans
        
        N1 = 0
        for m in range(0,N):
            if np.imag(roetter[m]) == 0:
                N1 = m
        if N1 < N:
            roetter[N1+1:N] = np.sort(roetter[N1+1:N])
        
        roetter = roetter-2*1j*np.imag(roetter)
        SERA = roetter
        
    #=========================================
    #   RESIDUE IDENTIFICATION
    #=========================================
    if opts.skip_res != 1:
        
        #We now calculate SER for f, using the modified zeros of sigma as new poles:
        
        LAMBD = roetter
        
        cindex = np.zeros([1,N])
        for m in range(0,N):
            if np.imag(LAMBD[m]) != 0:
                if m == 1:
                    cindex[m] = 1
                else:
                    if cindex[m-1] == 0 or cindex[m-1] == 2:
                        cindex[m] = 1
                        cindex[m+1] = 2
                    else:
                        cindex[m] = 2
        
        #====================================================================
        # We now calculate the SER for f (new fitting), using the above calc
        # zeros as known poles :
        #====================================================================
        
        if opts.asymp == 1:
            A = np.zeros([2*Ns,N])
            BB = np.zeros([2*Ns,Nc])
        elif opts.asymp == 2:
            A = np.zeros([2*Ns,N+1])
            BB = np.zeros([2*Ns,Nc])
        else:
            A = np.zeros([2*Ns,N+2])
            BB = np.zeros([2*Ns,Nc])
        
        Dk = np.zeros([Ns,N])
        for m in range(0,N):
            if cindex[m] == 0:  #Real pole
                Dk[:,m] = 1/(s-LAMBD[m])
            elif cindex[m] == 1:    #Complex pole, first part
                Dk[:,m] = weight/(s-LAMBD[m]) + weight/(s-np.transpose(LAMBD[m]))
                Dk[:,m+1] = 1j*weight/(s-LAMBD[m]) - 1j*weight/(s-np.transpose(LAMBD[m]))
        
        if common_weight == 1:
            
            Dk = np.zeros([Ns,N])
            for m in range(0,N):
                if cindex[m] == 0:  #Real pole
                    Dk[:,m] = weight/(s-LAMBD[m])
                elif cindex[m] == 1:    #Complex pole, first part
                    Dk[:,m] = weight/(s-LAMBD[m]) + weight/(s-np.transpose(LAMBD[m]))
                    Dk[:,m+1] = 1j*weight/(s-LAMBD[m]) - 1j*weight/(s-np.transpose(LAMBD[m]))
            
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
                BB[0:Ns,m] = weight*np.transpose(f[m,:])
            A[Ns+1:2*Ns,:] = np.imag(A[0:Ns,:])
            A[0:Ns,:] = np.real(A[0:Ns,:])
            BB[Ns+1:2*Ns,:] = np.imag(BB[0:Ns,:])
            BB[0:Ns,:] = np.real(BB[0:Ns,:])
            
            if opts.asymp == 2:
                A[0:Ns,N+1] = A[0:Ns,N+1]
            elif opts.asymp == 3:
                A[0:Ns,N+1] = A[0:Ns,N+1]
                A[Ns+1:2*Ns] = A[Ns+1:2*Ns,N+2]
                
            #clear Escale
            Escale = np.zeros([1,len(A[0,:])])
            for col in range(0,len(A[0,:])):
                Escale[col] = lin.norm(A[:,col],2)
                A[:,col] = A[:,col]/Escale[col]
            X = A/BB
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
        
        Dk = np.zeros([Ns,N])
        for m in range(0,N):
            Dk[:,m] = 1/(s-SERA[m])
        for n in range(0,Nc):
            fit[n,:] = np.transpose(Dk*np.transpose(SERC[n,:]))
            if opts.asymp == 2:
                fit[n,:] = fit[n,:] + SERD[n]
            elif opts.asymp == 3:
                fit[n,:] = fit[n,:] + SERD[n]+np.tranpose(s)*SERE[n]
        
        fit = np.transpose(fit)
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
        
        cindex = np.zeros([1,N])
        for m in range(0,N):
            if np.imag(A[m,m]) != 0:
                if m == 1:
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
    class SER():
        A = A
        B = B
        C = C
        D = D
        E = E
        
    return SER, poles, rmserr, fit, opts