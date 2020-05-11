import numpy as np
from numpy import log10, ceil, floor
from math import sqrt
from vectfit3 import vectfit3, DefaultVect3Opts, SERClass
import numpy.linalg as lin

class DefVFOpts():
    N = 10
    poletype = 'linlogcmplx'
    nu = 1e-3
    Niter1 = 7
    Niter2 = 0
    weight = []
    weightparam = 1
    asymp = 2
    stable = 1
    relaxed = 1
    plot = 1
    logx = 1
    logy = 1
    errplot = 1
    phaseplot = 1
    screen = 1
    cmplx_ss = 1
    remove_HFpoles = 0
    factor_HF = 1.1
    passive_DE = 0
    passive_DE_TOLD = 1e-6
    passive_DE_TOLE = 1e-16
    spy1 = 0
    spy2 = 0



def VFDriver(bigH, s, poles, opts):
    
    vect3opts = DefaultVect3Opts()
    
    Niter1 = opts.Niter1
    Niter2 = opts.Niter2
    weightparam = opts.weightparam
    remove_HFpoles = opts.remove_HFpoles
    factor_HF = opts.factor_HF
    passive_DE = opts.passive_DE
    passive_DE_TOLD = opts.passive_DE_TOLD
    passive_DE_TOLE = opts.passive_DE_TOLE
    
    fit1 = []
    fit2 = []
    fit3 = []
    Ns = len(s[0])
    
    if len(poles) == 0:
        if opts.N == 0:
            print('===> ERROR: You did not specify a value for opts.N (fitting order).')
            print('Solution: either specify value for opts.N, or provide initial poles in array poles.')
            return
        N = opts.N
        oldpoletype = opts.poletype
        if N < 6:
            if opts.poletype == 'linlogcmplx':
                opts.poletype = 'logcmplx'
        nu = opts.nu
        if len(opts.poletype) == 8: #N = 10, 5 complex conjugate pairs 
            if opts.poletype == 'lincmplx': #Complex, linearly spaced starting poles
                bet = np.linspace(s[0,0]/1j,s[0,-1]/1j,int(floor(N/2)))
                for n in range (0,len(bet)): #n (length of bet) = 5, N (ORDER) = 10
                    alf = -nu*bet[n]
                    piece1 = alf-1j*bet[n] #Complex conjugate
                    piece2 = alf+1j*bet[n] #Complex conjugate
                    poles = np.append(poles,piece1)
                    poles = np.append(poles,piece2)
                poles = poles.reshape(1,-1)
            elif opts.poletype == 'logcmplx': #Complex, log spaced starting poles
                bet = np.logspace(log10(s[0,0]/1j),log10(s[0,-1]/1j),int(floor(N/2)))
                poles = np.zeros([1,N], dtype = 'complex128')
                for n in range(0,len(bet)):
                    alf = -nu*bet[n]
                    piece1 = alf-1j*bet[n] #Complex conjugate
                    piece2 = alf+1j*bet[n] #Complex conjugate
                    poles = np.append(poles,piece1)
                    poles = np.append(poles,piece2)
                poles = poles.reshape(1,-1)
            else:
                print('ERROR: Illegal value for opts.poletype')
                print('Valid input: ''lincmplex'' and ''logcmplx''')
                print('Given input:')
                print(opts.poletype)
                return
            
        elif len(opts.poletype) == 11:
            if opts.poletype == 'linlogcmplx':
                poles1 = []
                bet = np.linspace(s[0,0]/1j, s[0,-1]/1j, int(ceil((N-1)/4)))
                for n in range(0,len(bet)):
                    alf = -nu*bet[n]
                    piece1 = alf-1j*bet[n] #Complex conjugate
                    piece2 = alf+1j*bet[n] #Complex conjugate
                    poles1 = np.append(poles1,piece1)
                    poles1 = np.append(poles1,piece2)
                bet = np.logspace(log10(s[0,0]/1j), log10(s[0,-1]/1j), int(2+floor(N/4)))
                bet = bet[1:-1]
                poles2 = []
                for n in range(0,len(bet)):
                    alf = -nu*bet[n]
                    piece1 = alf-1j*bet[n] #Complex conjugate
                    piece2 = alf+1j*bet[n] #Complex conjugate
                    poles2 = np.append(poles2,piece1)
                    poles2 = np.append(poles2,piece2)
                poles = np.append(poles1, poles2)
                poles = poles.reshape(1,-1)
                
            else:
                print('ERROR: Illegal value for opts.poletype')
                print('Valid input: ''lincmplex'' and ''logcmplx''')
                print('Given input:')
                print(opts.poletype)
                return
        
        if len(poles[0]) < N:
            if opts.poletype == 'lincmplx':
                pole_extra = -(s[0,0]/1j+s[0,-1]/1j)/2
            elif opts.poletype == 'logcmplx' or opts.poletype == 'linlogcmplx':
                pole_extra = -10**((log10(s[0,0]/1j)+log10(s[0,-1]/1j))/2)
            poles = np.append(poles, pole_extra)
        opts.poletype = oldpoletype
    
    
    Nc = len(bigH[:,0,0])
    Ns = len(s[0])
    
    f = np.zeros([sum(range(Nc+1)),len(bigH[0,0,:])], dtype = 'complex128')
    
    if opts.screen == 1:
        print('==================== START ====================')
    tell = -1
    for col in range(0,Nc):
        for row in range(col,Nc):
            tell = tell + 1
            f[tell,:] = np.transpose(np.squeeze(bigH[row,col,:])) #Stacking elements into a single vector
        
    nnn = tell
    
    #Fitting options
    vect3opts.spy1 = 0
    vect3opts.spy2 = 0
    vect3opts.skip_pole = 0
    vect3opts.skip_res = 1
    vect3opts.legend = 1
    
    oldspy2 = vect3opts.spy2
    opts.spy2 = 0
    if Nc == 1:
        f_sum = f
    if Nc > 1:  #For the multi terminal case
        #Forming columns sum and associated LS weight
        f_sum = np.zeros([1,f.shape[1]])
        tell = -1
        for row in range(0,Nc):
            for col in range(row,Nc):
                tell = tell + 1
                if weightparam == 1 or weightparam == 4 or weightparam == 5:
                    f_sum = f_sum + f[tell,:]
                elif weightparam == 2:
                    f_sum = f_sum + f[tell,:]/lin.norm(f[tell,:])
                elif weightparam == 3:
                    f_sum = f_sum + f[tell,:]/sqrt(lin.norm(f[n,:]))
    
    
    #Creating LS weight
    if len(opts.weight) == 0: #Automatic specification of weight
        if weightparam == 1: #1 for all elements in LS problem
            weight = np.ones([1,Ns])
            weight_sum = np.ones([1,Ns])
        elif weightparam == 2: #Individual element weighting
            weight = 1/abs(f)
            weight_sum = 1/abs(f_sum)
        elif weightparam == 3: #Individual element weighting
            weight = 1/sqrt(abs(f))
            weight_sum = 1/sqrt(abs(f_sum))
        elif weightparam == 4: #Common weighting for all matrix elements
            for k in range(0,Ns):
                weight[k] = 1/lin.norm(f[:,k])
            weight_sum = weight
        elif weightparam == 5: #Common weighting for all mamtrix elements
            for k in range(0,Ns):
                weight[k] = 1/sqrt(lin.norm(f[:,k]))
            weight_sum = weight
        else:
            print('ERROR in mtrxVectfit: Illegal value for opts.weight')
            return
    else:
        weight = np.zeros([nnn,Ns])
        tell = -1
        for row in range(0,Nc):
            for col in range(row,Nc):
                tell = tell + 1
                weight[tell,:] = np.squeeze(opts.weight[row,col,:])
        weight_sum = np.ones(Ns)
    
    if Nc > 1:   #Will do only for multi-terminal case
        if opts.screen == 1:
            print('****Calculating improved initial poles by fitting column sum ...')
        for itr in range(0,Niter1):
            if opts.screen == 1:
                print('Iter ', str(itr))
            [SER,poles,rmserr,fit] = vectfit3(f_sum,s,poles,weight_sum,vect3opts)
    if opts.screen == 1:
        print('****Fitting column ...')
    vect3opts.skip_res = 1
    for itr in range(0,Niter2):
        if opts.screen == 1:
            print('Iter ', str(itr))
        if itr == Niter2:
            vect3opts.skip_res = 0
        [SER,poles,rmserr,fit1] = vectfit3(f,s,poles,weight,vect3opts)
    if Niter2 == 0:
        vect3opts.skip_res = 0
        vect3opts.skip_pole = 1
        [SER,poles,rmserr,fit1] = vectfit3(f,s,poles,weight,vect3opts) #HAVE STOPPED HERE
    
    #========================================
    #   Throwing out high frequency poles
    fit2 = fit1
    if remove_HFpoles == 1:
        if opts.screen == 1:
            print('****Throwing out high frequency poles: ...')
        ind = np.where(abs(poles)>factor_HF*abs(s[-1]))
        poles[ind] = []
        N = len(poles)
        if opts.screen == 1:
            print('****Refitting residues: ...')
        DefVFOpts.skip_pole = 1
        [SER,poles,rmserr,fit2] = vectfit3(fit1,s,poles,weight,DefVFOpts)
    #=========================================
    
    #=========================================
    if passive_DE == 1 and DefVFOpts.asymp > 1:
        if opts.screen == 1:
            if DefVFOpts.asymp == 2:
                print('****Enforcing positive realness for D...')
            elif DefVFOpts.asymp == 3:
                print('****Enforcing positive realness for D, E...')
        tell = 0
        DD = np.zeros(Nc)
        EE = np.zeros(Nc)
        for col in range(0,Nc):
            for row in range(col,Nc):
                tell = tell + 1
                DD[row,col] = SER.D[tell]
                EE[row,col] = SER.E[tell]
        #Calculating Dmod, Emod:
        [V,L] = lin.eig(DD)
        for n in range(0,Nc):
            if L[n,n] < 0:
                L[n,n] = passive_DE_TOLD
        DD = V*L*V**-1
        [V,L] = lin.eig(EE)
        for n in range(0,Nc):
            if L[n,n] < 0:
                L[n,n] = passive_DE_TOLE
        EE = V*L*V**-1
        tell = 0
        #Calculating fmod
        Emod = np.zeros(Nc)
        Dmod = np.zeros(Nc)
        fmod = np.zeros(Nc)
        for col in range (0,Nc):
            for row in range(col,Nc):
                tell = tell + 1
                Dmod[tell] = DD[row,col]
                Emod[tell] = EE[row,col]
                fmod[tell,:] = fit2[tell,:] - Dmod[tell] - s*Emod[tell]
        if opts.screen == 1:
            if vect3opts.asymp == 2:
                print('****Refitting C while enforcing D=0...')
            elif vect3opts.asymp == 3:
                print('****Refitting C while enforcing D=0, E=0...')
        vect3opts.skip_pole = 0
        vect3opts.asymp = 1
        for itr in range(0,0):
            [SER,poles,rmserr,fit3]=vectfit3(fmod,s,poles,weight,vect3opts)
        SER.D = Dmod
        SER.E = Emod
        for tell in range(0,len(fit3[:,0])):
            fit3[tell,:] = fit3[tell,:] + SER.D[tell] + s*SER.E[tell]
    
    if Nc > 1:
        if opts.screen == 1:
            print('****Transforming model of lower matrix triangle into state-space model of full matrix...')
        SER.A = SER.A*np.eye(N)
        SER = tri2full(SER) #ERRORS: A[31,30], B[31,3], 
    
    if opts.screen == 1:
        print('****Generating pole-residue model...')
    [R,a] = ss2pr(SER.A,SER.B,SER.C)
    SER.R = R
    SER.poles = a
    
    #rmserror of fitting:
    if len(fit3) != 0:
        fit = fit3
    elif len(fit2) != 0:
        fit = fit2
    elif len(fit1) != 0:
        fit = fit1
    fit = np.transpose(fit)
    diff = fit - f
    rmserr = sqrt(sum(sum(abs(diff**2))))/sqrt(nnn*Ns)
    
    vect3opts.spy2 = oldspy2
    if vect3opts.spy2 == 1:
        if opts.screen == 1:
            print('****Plotting of results')
        
        #freq = s/(s*pi*1j)
            
        ###PLOTTING, REVISIT
    
    if opts.screen == 1:
        print('===================== END =====================')
    
    bigHfit = np.zeros([Nc,Nc,Ns], dtype = 'complex128')
    tell = -1
    for row in range(0,Nc):
        for col in range(row,Nc):
            tell = tell + 1
            bigHfit[row,col,:] = fit[tell,:]
            if row != col:
                bigHfit[col,row,:] = fit[tell,:]
    
    return SER, rmserr, bigHfit, opts
    

#==================================================

def tri2full(SER):
    
    import scipy
    
    A = SER.A
    B = SER.B
    C = SER.C
    D = SER.D
    E = SER.E
    
    tell = 0
    for k in range(1,10000): #1e4
        tell = tell + k
        if tell == len(D):
            Nc = k
            break
    N = len(A)
    tell = -1
    CC = np.zeros([Nc,Nc*N], dtype = 'complex128')
    AA = None
    BB = None
    DD = np.zeros([Nc,Nc], dtype = 'complex128')
    EE = np.zeros([Nc,Nc], dtype = 'complex128')
    for col in range(0,Nc):
        if col == 0:
            AA = A
            BB = B
        else:
            AA = scipy.linalg.block_diag(AA,A)
            BB = scipy.linalg.block_diag(BB,B)
        for row in range(col,Nc):
            tell = tell + 1
            DD[row,col] = D[tell]
            EE[row,col] = E[tell]
            tmpci1 = (col)*N
            tmpci12 = (col+1)*N
            tmpci2 = (row)*N
            tmpci22 = (row+1)*N
            CC[row,tmpci1:tmpci12] = C[tell,:]
            CC[col,tmpci2:tmpci22] = C[tell,:]
            #CC[row,(col-1)*N+1:col*N] = C[tell,:]
            #CC[col,(row-1)*N+1:row*N] = C[tell,:]
    DD = DD + np.transpose(DD-np.diag(np.diag(DD)))
    EE = EE + np.transpose(EE-np.diag(np.diag(DD)))
    
    SER = SERClass(AA,BB,CC,DD,EE)

    
    return SER


#============================================================

def ss2pr(A,B,C):
    
    #Converting real-only state-space model into complex model
    tmp1 = np.diag(A)
    tmp2 = np.diag(tmp1)
    tmp3 = A - tmp2
    tmp4 = abs(tmp3)
    tmp5 = np.max(tmp4)
    tmp6 = np.max(tmp5)
    if tmp6 != 0:
        #if max(max(abs(A-np.diag(np.diag(A))))) != 0:
        #errflag = 0
        for m in range(0,len(A)-1):
            if A[m,m+1] != 0:
                A[m,m] = A[m,m] +1j*A[m,m+1]
                A[m+1,m+1] = A[m+1,m+1] - 1j*A[m,m+1]
                
                B[m,:] = (B[m,:]+B[m+1,:])/2
                B[m+1,:] = B[m,:]
                
                C[:,m] = C[:,m] + 1j*C[:,m+1]
                C[:,m+1] = np.conj(C[:,m])
                
    #Converting complex state-space model into pole-residue model
    Nc = len(C[:,0])
    N = int(len(A)/Nc)
    R = np.zeros([Nc,Nc,N], dtype = 'complex128')
    for m in range(0,N):
        Rdum = np.zeros(Nc, dtype = 'complex128')
        for n in range(0,Nc):
            ind = (n)*N+m
            Rdum = Rdum + C[:,ind]*B[ind,:]
        R[:,:,m] = Rdum
    a = np.diag(A[0:N,0:N])
    
    return R, a
    
