def RPDriver(SER, s, opts):
    MPopts.auxflag = 1
    MPopts.solver = 'QUADPROG'
    
    [SER] = pr2ss(SER)  #Convert model from pole-residue to state-space
    
    print('================== START ==================')
    
    class opts():
        parametertype = 'Y'
        Niter_out = 10
        Niter_in = 0
        TOLGD = 1e-6
        TOLE = 1e-12
        cmplx_ss = 1
        weightfactor = 0.001
        weightparam = 1
        method = 'FRP'
        colinterch = 1
        outputlevel = 1
        weight = []
        
    if opts.method == 'FMP' and opts.parametertype == 'S':
        print('ERROR in RPDriver: FMP cannot be used together with S-parameters')
        return
    
    colinterch = opts.colinterch
    MPopts.TOLGD = opts.TOLGD
    MPopts.TOLE = opts.TOLE
    MPopts.weightfactor = opts.weightfactor
    MPopts.weight = opts.weight
    MPopts.outputlevel = opts.outputlevel
    
    #Check if line 163 is necessary
    
    if opts.parametertype == 'Y':
        print('*** Y-PARAMETERS ***')
    elif opts.parametertype == 'S':
        print('*** S-PARAMETERS ***')
    
    #Check lines 174-182
    break_outer = 0
    olds3 = []
    
    SER0 = SER
    Nc = len(SER.D)
    
    Niter_out = opts.Niter_out
    Niter_in = opts.Niter_in
    
    #============================================
    #   Plotting eigenvalues of orinigal model
    #============================================
    
    #Check lines 195-236
    
    outputlevel = opts.outputlevel
    t = [0,0,0,0]
    
    #============================================
    #   Passivity enforcement
    #============================================
    
    QP.first = 1
    QPopts = []
    SER1 = SER0
    
    for iter_out in range(0,Niter_out):
        if break_outer == 1:
            SER0 = SER1
            break
        s3 = []
        
        for iter_in in range(0,0,Niter_in+1):
            s2 = []
            SERflag = 1
            if outputlevel == 1:
                
    
    
    