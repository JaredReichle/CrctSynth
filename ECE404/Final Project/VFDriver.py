def VFDriver(bigH, s, poles, opts):

    import numpy as np
    from numpy import log10, ceil, floor
    from math import sqrt, pi
    import vectfit3
    
    class defa():
        N = []
        poletype = 'lincmplx'
        nu = 1e-3
        Niter1 = 4
        Niter2 = 4
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
    
    A = fieldnames(defa)
    for m in range(0,len(A)):
        if isfield(opts,A[m]) != bool(true):
            dum = char(A[m])
            dum2 = getfield(defa,dum)
            opts = setfield(opts,dum,dum2)
    
    VF.asymp = opts.asymp
    if opts.stable == 1:
        VF.stable = 1
    else:
        VF.stable = 0
    if opts.relaxed == 1:
        VF.relax = 1
    else:
        VF.relax = 0
    if opts.plot == 1:
        VF.spy2 = 1
    else:
        VF.spy2 = 0
    if opts.logx == 1:
        VF.logx = 1
    else:
        VF.logx = 0
    if opts.logy == 1:
        VF.logy = 1
    else:
        VF.logy = 0
    if opts.errplot == 1:
        VF.errplot = 1
    else:
        VF.errplot = 0
    if opts.phaseplot == 1:
        VF.phaseplot = 1
    else:
        VF.phaseplot = 0
    if opts.cmplx_ss == 1:
        VF.complx_ss = 1
    else:
        VF.complx_ss = 0
        
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
    Ns = len(s)
    
    if isempty(poles):
        if isempty(opts.N):
            print('===> ERROR: You did not specify a value for opts.N (fitting order).')
            print('Solution: either specify value for opts.N, or provide initial poles in array poles.')
            return
        N = opts.N
        oldpoletype = opts.poletype
        if N < 6:
            if strcmp(opts.poletype, 'linlogcmplx'):
                opts.poletype = 'logcmplx'
        nu = opts.nu
        if len(opts.poletype) == 8:
            if opts.poletype == 'lincmplx': #Complex, linearly spaced starting poles
                bet = np.linspace(s[0]/1j,s[Ns]/1j,floor(N/2))
                poles = []
                for n in range (0,len(bet)):
                    alf = -nu*bet[n]
                    poles = [poles, (alf-1j*bet[n]), (alf+1j*bet[n])]
            elif opts.poletype == 'logcmplx': #Complex, log spaced starting poles
                bet = np.logspace(log10(s[0]/1j),log10(s[Ns]/1j),floor(N/2))
                poles = []
                for n in range(0,len(bet)):
                    alf = -nu*bet[n]
                    poles = [poles, (alf-1j*bet[n]), (alf+1j*bet[n])]
            else:
                print('ERROR: Illegal value for opts.poletype')
                print('Valid input: ''lincmplex'' and ''logcmplx''')
                print('Given input:')
                print(opts.poletype)
                return
            
        elif len(opts.poletype) == 11:
            if opts.poletype == 'linlogcmplx':
                bet = np.linspace(s[0]/1j, s[Ns]/1j, ceil((N-1)/4))
                poles1 = []
                for n in range(0,len(bet)):
                    alf = -nu*bet[n]
                    poles1 = [poles1, (alf-1j*bet[n]), (alf+1j*bet[n])]
                bet = np.logspace(log10(s[0]/1j), log10(s[Ns]/1j), 2+floor(N/4))
                bet[0] = []
                bet[-1] = []
                poles2 = []
                for n in range(0,len(bet)):
                    alf = -nu*bet[n]
                    poles2 = [poles2, (alf-1j*bet[n]), (alf+1j*bet[n])]
                poles = [poles1, poles2]
                
            else:
                print('ERROR: Illegal value for opts.poletype')
                print('Valid input: ''lincmplex'' and ''logcmplx''')
                print('Given input:')
                print(opts.poletype)
                return
        
        if len(poles) < N:
            if strcmp(opts.poletype, 'lincmplx'):
                pole_extra = -(s[0]/1j+s[-1]/1j)/2
            elif strcmp(opts.poletype, 'logcmplx') or strcmp(opts.poletype,'linlogcmplx'):
                pole_extra = -10**((log10(s[0]/1j)+log10(s[-1]/1j))/2)
            poles = [poles, pole_extra]
        
        opts.poletype = oldpoletype
        
    Nc = len(bigH[:,0,0])
    Ns = len(s)
    
    if opts.screen == 1:
        print('==================== START ====================')
    tell = 0
    for col in range(0,Nc):
        for row in range(col,Nc):
            tell = tell + 1
            f[tell,:] = np.transpose(np.squeeze(bigH[row,col,:])) #Stacking elements into a single vector
        
    nnn = tell
    
    #Fitting options
    VF.spy1 = 0
    VF.skip_pole = 0
    VF.skip_res = 1
    VF.legend = 1
    
    oldspy2 = VF.spy2
    VF.spy2 = 0
    if Nc == 1:
        f_sum = f
    if Nc > 1:  #For the multi terminal case
        #Forming columns sum and associated LS weight
        f_sum = 0
        tell = 0
        for row in range(0,Nc):
            for col in range(row,Nc):
                tell = tell + 1
                if weightparam == 1 or weightparam == 4 or weightparam == 5:
                    f_sum = f_sum+f[tell,:]
                elif weightparam == 2:
                    f_sum = f_sum + f[tell,:]/norm(f[tell,:])
                elif weightparam == 3:
                    f_sum = f_sum + f[tell,:]/sqrt(norm(f[n,:]))
    
    #Creating LS weight
    if isempty(opts.weight): #Automatic specification of weight
        if weightparam == 1: #1 for all elements in LS problem
            weight = np.ones(Ns)
            weight_sum = np.ones(Ns)
        elif weightparam == 2: #Individual element weighting
            weight = 1/abs(f)
            weight_sum = 1/abs(f_sum)
        elif weightparam == 3: #Individual element weighting
            weight = 1/sqrt(abs(f))
            weight_sum = 1/sqrt(abs(f_sum))
        elif weightparam == 4: #Common weighting for all matrix elements
            for k in range(0,Ns):
                weight[k] = 1/norm(f[:,k])
            weight_sum = weight
        elif weightparam == 5: #Common weighting for all mamtrix elements
            for k in range(0,Ns):
                weight[k] = 1/sqrt(norm(f[:,k]))
            weight_sum = weight
        else:
            print('ERROR in mtrxVectfit: Illegal value for opts.weight')
            return
    else:
        weight = np.zeros([nnn,Ns])
        tell = 0
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
                print('Iter ', num2str(itr))
            [SER,poles,rmserr,fit] = vectfit3(f_sum,s,poles,weight_sum,VF)
    if opts.screen == 1:
        print('****Fitting column ...')
    VF.skip_res = 1
    for itr in range(0,Niter2):
        if opts.screen == 1:
            print('Iter ', num2str(itr))
        if itr == Niter2:
            VF.skip_res = 0
        [SER,poles,rmserr,fit1] = vectfit3(f,s,poles,weight,VF)
    if Niter2 == 0:
        VF.skip_res = 0
        VF.skip_pole = 1
        [SER,poles,rmserr,fit1] = vectfit3(f,s,poles,weight,VF)
    
    #========================================
    #   Throwing out high frequency poles
    fit2 = fit1
    if remove_HFpoles == 1:
        if opts.screen == 1:
            print('****Throwing out high frequency poles: ...')
        ind = find(abs(poles)>factor_HF*abs(s[-1]))
        poles[ind] = []
        N = len(poles)
        if opts.screen == 1:
            print('****Refitting residues: ...')
        VF.skip_pole = 1
        [SER,poles,rmserr,fit2] = vectfit3(fit1,s,poles,weight,VF)
    #=========================================
    
    #=========================================
    if passive_DE == 1 and VF.asymp > 1:
        if opts.screen == 1:
            if VF.asymp == 2:
                print('****Enforcing positive realness for D...')
            elif VF.asymp == 3:
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
        [V,L] = la.eig(DD)
        for n in range(0,Nc):
            if L[n,n] < 0:
                L[n,n] = passive_DE_TOLD
        DD = V*L*V**-1
        [V,L] = la.eig(EE)
        for n in range(0,Nc):
            if L[n,n] < 0:
                L[n,n] = passive_DE_TOLE
        EE = V*L*V**-1
        tell = 0
        #Calculating fmod
        Emod = np.zeros(Nc)
        Dmod = np.zeros(Nc)
        for col in range (0,Nc):
            for row in range(col,Nc):
                tell = tell + 1
                Dmod[tell] = DD[row,col]
                Emod[tell] = EE[row,col]
                fmod[tell,:] = fit2[tell,:] - Dmod[tell] - s*Emod[tell]
        if opts.screen == 1:
            if VF.asymp == 2:
                print('****Refitting C while enforcing D=0...')
            elif VF.asymp == 3:
                print('****Refitting C while enforcing D=0, E=0...')
        VF.skip_pole = 0
        VF.asymp = 1
        for itr in range(0,0):
            [SER,poles,rmserr,fit3]=vectfit3(fmod,s,poles,weight,VF)
        SER.D = Dmod
        SER.E = Emod
        for tell in range(0,len(fit3[:,0])):
            fit3[tell,:] = fit3[tell,:] + SER.D[tell] + s*SER.E[tell]
    
    if Nc > 1:
        if opts.screen == 1:
            print('****Transforming model of lower matrix triangle into state-space model of full matrix...')
        [SER] = tri2full(SER) #FIND PYTHON ROUTINE
    
    if opts.screen == 1:
        print('****Generating pole-residue model...')
    [R,a] == ss2pr(SER.A,SER.B,SER.C) #FIND PYTHON ROUTINE
    SER.R = R
    SER.poles = a
    
    #rmserror of fitting:
    if isempty(fit3) == 0:
        fit = fit3
    elif isempty(fit2) == 0:
        fit = fit2
    elif isempty(fit1) == 0:
        fit = fit1
    diff = fit - f
    rmserr = sqrt(sum(sum(abs(diff**2))))/sqrt(nnn*Ns)
    
    VF.spy2 = oldspy2
    if VF.spy2 == 1:
        if opts.screen == 1:
            print('****Plotting of results')
        
        freq = s/(s*pi*1j)
            
        ###PLOTTING, REVISIT
    
    if opts.screen == 1:
        print('===================== END =====================')
    
    bigHfit = np.zeros([Nc,Nc,Ns])
    tell = -1
    for row in range(0,Nc):
        for col in range(row,Nc):
            tell = tell + 1
            bighfit[row,col,:] = fit[tell,:]
            if row != col:
                bigHfit[col,row,:] = fit[tell,:]
    
    return SER, rmserr, bigHfit, opts
    
def tri2full(SER):
    
    import scipy.linalg as lin
    
    A = SER.A
    B = SER.B
    C = SER.C
    D = SER.D
    E = SER.E
    
    tell = 0
    for k in range(0,1e4):
        tell = tell + k
        if tell == len(D):
            Nc = k
            break
    N = len(A)
    tell = -1
    CC = np.zeros([Nc,Nc*N])
    AA = []
    BB = []
    for col in range(0,Nc):
        AA = lin.block_diag(AA,A)
        BB = lin.block_diag(BB,B)
        for row in range(col,Nc):
            tell = tell + 1
            DD[row,col] = D[tell]
            EE[row,col] = E[tell]
            CC[row,(col-1)*N+1:col*N] = C[tell,:]
            CC[col,(row-1)*N+1:row*N] = C[tell,:]
    DD = DD + np.transpose(DD-lin.diag(lin.diag(DD)))
    EE = EE + np.transpose(EE-lin.diag(lin.diag(DD)))
    
    SER2.A = AA
    SER2.B = BB
    SER2.C = CC
    SER2.D = DD
    SER2.E = EE
    
    return SER2

def ss2pr(A,B,C):
    
    import numpy as np
    import scipy.linalg as lin
    
    #Converting real-only state-space model into complex model
    if max(max(abs(A-lin.diag(lin.diag(A))))) != 0:
        errflag = 0
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
    N = len(A)/Nc
    R = np.zeros([Nc,Nc,N])
    for m in range(0,N):
        Rdum = np.zeros(Nc)
        for n in range(0,Nc):
            ind = (n-1)*N+m
            Rdum = Rdum + C[:,ind]*B[ind,:]
        R[:,:,m] = Rdum
    a = np.full(lin.diag(A[0:N,0:N]))
    
    return R, a
    