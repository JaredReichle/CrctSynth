import numpy as np
import numpy.linalg as lin
from math import pi
import scipy.linalg


#from pr2ss import pr2ss
from violextremaY import violextremaY
from FRPY import FRPY, FRPOpts
from intercheig import intercheig
from rot import rot
from fitcalcABCDE import fitcalcABCDE

class DefRPOpts():
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
    s_pass = np.transpose(2*pi*1j*np.linspace(0,2E5,1001))
    ylim = [-2e-3, 2e-3]


#LINE 256 MATLAB    
class QP():
    first = 1

def RPDriver(SER, s, opts):
 
    #SER = pr2ss(SER)  #Convert model from pole-residue to state-space

    print('================== START ==================')
    
    frpopts = FRPOpts()
    
    frpopts.TOLGD = opts.TOLGD
    frpopts.TOLE = opts.TOLE
    frpopts.weightfactor = opts.weightfactor
    frpopts.weight = opts.weight
    frpopts.outputlevel = opts.outputlevel
    
    colinterch = opts.colinterch
    
    print('*** Y-PARAMETERS ***')
    
    #Check lines 174-182 for plotting
    
    break_outer = 0
    
    #olds3 = []
    
    SER0 = SER
    Nc = len(SER.D)
    
    Niter_out = opts.Niter_out
    Niter_in = opts.Niter_in
    
    #============================================
    #   Plotting eigenvalues of orinigal model
    #============================================
    
    #Check lines 195-236 - Revisit
    
    outputlevel = opts.outputlevel
    
    #============================================
    #   Passivity enforcement
    #============================================
    
    QP.first = 1
    #QPopts = []
    
    SER1 = SER0
    
    for iter_out in range(0,Niter_out):
        if break_outer == 1:
            SER0 = SER1
            break
        s3 = []
        
        for iter_in in range(0,Niter_in+1):
            s2 = []
            SERflag = 0
            if outputlevel == 1:
                print('  ')
                print(str(iter_out), '  ', str(iter_in), '  ')
            
            if iter_in == 0:
                wintervals = pass_check_Y(SERflag, SER.A,SER.B,SER.C,SER.D,colinterch)
                #[wintervals] = pass_check_Y(SERflag, SER.poles,[],SER1.R,SER1.D,colinterch)
                
                #Still have violating intervals when none are present
                
                if len(wintervals) > 0:
                    if outputlevel == 1:
                        print('N.o. violating intervals: ', str(len(wintervals[0,:])))
                
                if len(wintervals) == 0 and all(lin.eig(SER1.D) >= 0) and all(lin.eig(SER1.E)):
                    SER0 = SER1
                    break_outer = 1
                    break
                
                [s_viol,g_pass,ss] = violextremaY(SERflag,np.transpose(wintervals),SER.poles,[],SER1.R, SER1.D,colinterch)
                
                s2 = [ np.transpose(s_viol)]
                s2 = np.sort(s2)
                if len(s2) == 0 and all(lin.eig(SER1.D)>0):
                    break
            
            if iter_in == 0:
                if outputlevel == 1:
                    if min(g_pass) < 0:
                        print('Max violation, eig(G) : ', str(g_pass), ' @ ', str((ss)/(2*pi*1j)))
                    else:
                        print('Max violation, eig(G) : None')
                    
                    if min(lin.eig(SER0.D)) < 0:
                        print('Max violation, eig(G) : ', str(min(lin.eig(SER1.D))))
                    else:
                        print('Max violation, eig(G) : None')
                    
                    if min (lin.eig(SER0.E)) < 0:
                        print('Max violation, eig(G) : ', str(min(lin.eig(SER1.E))))
                    else:
                        print('Max violation, eig(G) ; None')
                if outputlevel != 1:
                    min1 = min(g_pass)
                    min2 = min(lin.eig(SER1.D))
                    print('Max violation : ', str(min([min1,min2])))
                    if min(lin.eig(SER0.E)) < 0:
                        print('Max violation, E : ',str(min(lin.eig(SER1.E))))
                    
            
            if outputlevel == 1:
                print('Passivity Enforcement...')
            #if opts.method == 'FMP':
                #[SER1, MPopts] = FMP(SER0,s,s2,s3,MPopts)
                #No such FMP routine
                if opts.method == 'FRP': #Can't find FMP routine
                    [SER1, MPopts] = FRPY(SER0,s,s2,s3,MPopts)
                else:
                    print('****** ERROR #1 in RPDriver.py')
            
            #if plotte == 1:
                
                #PLOTTING FUNCTIONS
                
            #end plotte
            #Confusing. Check below
            if iter_in != Niter_in + 1:
                [wintervals] = pass_check_Y(SERflag,SER1.poles,[],SER1.R,SER1.D)
                [s_viol] = violextremaY(SERflag,np.transpose(wintervals),SER1.poles,[],SER1.R,SER1.D,colinterch)
                
            #olds3 = s3
            s3 = [s3,s2,np.transpose(s_viol)]
            
            #Confusing. Check below
            if iter_in == Niter_in + 1:
                s3 = []
                s2 = []
                if plotte == 1:
                    EE0 = EE1
                SER0 = SER1
    
    #===========================================================
    #   Plotting eigenvalues of modified model (SERC1, SERD1)
    #===========================================================

    s_pass = opts.plot.s_pass

    EE1 = []
    if plotte == 1:
        oldT0 = []
        #tell = -1
        for k in range(0,len(s_pass)):
            Y = SER1.C*lin.diag((s_pass[k]*lin.I - lin.diag(SER1.A))**(-1))*SER1.B+SER1.D+s_pass[k]*SER1.E
            G = np.real(Y)
        [T0,D] = lin.eig(G)
        T0 = rot(T0)    #Find routine - minimizing phase angle of eigenvectors
        [T0,D] = intercheig(T0, oldT0, D, Nc, k) #Find routine
        oldT0 = T0
        EE1[:,k] = lin.diag(D)
        
        #Plotting items
    
    if wintervals.size == 0:
        if outputlevel == 1:
            print('   ')
        print('---> Passivity was successfully enforced')
        if outputlevel == 1:
            print('Max violation, eig(G) : None')
            print('Max violation, eig(D) : None')
            print('Max violation, eig(E) : None')
    else:
        print('Max violation, eig(G) : ', str(min(g_pass)))
        print('Max violation, eig(D) : ', str(min(lin.eig(SER0.D))))
        print('Max violation, eig(E) : ', str(min(lin.eig(SER0.E))))
        print('--> Iterations terminated before completing passivity enforcement.')
        print('Increase parameter opts.Niter_out')
    
    #Producing plot
    Ns = len(s)
    bigYfit = np.zeros([Nc,Nc,Ns])
    I = sci.sparse(np.ones([len(SER.A[:,0]),1]))
    for k in range(0,Ns):
        Y = SER1.C*lin.diag((s_pass[k]*I - lin.diag(SER1.A))**(-1))*SER1.B+SER1.D+s_pass[k]*SER1.E
        bigYfit[:,:,k] = Y
    
    #553
       
    if opts.cmplx_ss == 0:
        N = len(SER1.A)
        cindex = np.zeros([1,N])
        for m in range(0,N):
            if np.imag(SER1.A[m,m]) != 0:
                if m == 1:
                    cindex[m] = 1
                else:
                    if cindex[m-1] == 0 or cindex[m-1] == 2:
                        cindex[m] = 1
                        cindex[m+1] = 2
                    else:
                        cindex[m] = 2
            
        n = -1
        Ablock = []
        for m in range(0,N):
            n = n + 1
            if cindex[m] == 1:
                a = SER1.A[n,n]
                a1 = np.real(a)
                a2 = np.imag(a)
                c = SER1.C[:,n]
                c1 = np.real(c)
                c2 = np.imag(c)
                b = SER1.B[n,:]
                b1 = 2*np.real(b)
                b2 = -2*np.imag(b)
                Ablock[[a1,a2],[-a2,a1]]
                SER1.A[n:n+1,n:n+1] = Ablock
                SER1.C[:,n] = c1
                SER1.C[:,n+1] = c2
                SER1.B[n,:] = b1
                SER1.B[n+1,:] = b2
    
    print('================ END ================')
    
    return SER1, bigYfit, opts
    
def pass_check_Y(SERflag,A,B,C,D,colinterch):
    wintervals = []
    
    if SERflag == 1:
        Nc = len(D)
        N = len(A)
        tell = -1
        A = np.diag(A)
        CC = np.zeros([Nc,Nc*N], dtype = 'complex128')
        B = np.ones([N,1], dtype = 'complex128')
        for col in range(0,Nc):
            if col == 0:
                AA = A
                BB = B
            else:
                AA = scipy.linalg.block_diag(AA,A)
                BB = scipy.linalg.block_diag(BB,B)
            for row in range(col,Nc):               
                tmpci1 = (col)*N
                tmpci12 = (col+1)*N
                tmpci2 = (row)*N
                tmpci22 = (row+1)*N
                CC[row,tmpci1:tmpci12] = C[row,col,:]
                CC[col,tmpci2:tmpci22] = C[row,col,:]
                #CC[row,(col-1)*N+1:col*N] = C[row,col,:]
                #CC[col,(row-1)*N+1:row*N] = C[row,col,:]
        A = AA
        B = BB
        C = CC
        
    Acmplx = A
    Bcmplx = B
    Ccmplx = C
    Dcmplx = D
    
    if sum(sum(A-A)) == 0:
        N = len(A)
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
                b2 = -2*np.imag(b)
                Ablock = np.vstack(([a1,a2],[-a2,a1]))
                A[n:n+2,n:n+2] = Ablock
                C[:,n] = c1
                C[:,n+1] = c2
                B[n,:] = b1
                B[n+1,:] = b2
                
    N = len(A)
    Nc = len(D)
    tell = -1
    E=np.zeros([Nc,Nc])
    tmp1 = lin.eig(D)[0]
    tmp2 = tmp1 == 0
    tmp3 = any(tmp2)
    tmp4 = tmp3 != False
    if tmp4:
    #if sum(lin.eig(D) == 0) > 0:
        Ahat = np.linalg.solve(A,np.eye(N))
        Bhat = np.matmul(Ahat,B)
        Chat = np.matmul(C,Ahat)
        Ddum = np.matmul(Ahat,B)
        Dhat = D-np.matmul(C,Ddum)
        A = Ahat
        B = Bhat
        C = Chat
        D = Dhat
    
    tmp1 = D**(-1)
    tmp2 = np.matmul(B,tmp1)
    tmp3 = np.matmul(tmp2,C)
    tmp4 = tmp3 - A
    tmp5 = np.matmul(A,tmp4)
    
    S1 = tmp5
    
    wS1 = lin.eig(S1)[0]
    wS1 = np.sqrt(wS1)
    tmp1 = lin.eig(Dcmplx)[0]
    tmp2 = tmp1 == 0
    tmp3 = any(tmp2)
    tmp4 = tmp3 != False
    if tmp4:
        #if sum(lin.eig(Dcmplx)==0)>0:
        wS1 = 1/wS1
    ind = np.where(np.imag(wS1) == 0)
    wS1 = wS1[ind]
    sing_w = np.sort(wS1)
    if len(sing_w) == 0:
        sing_w = []
        intervals = []
        return
    A = Acmplx
    B = Bcmplx
    C = Ccmplx
    D = Dcmplx
    
    viol = np.zeros(len(sing_w))
    
    midw = np.zeros([len(sing_w),1])
    midw[0] = sing_w[0]/2
    midw[-1] = 2*sing_w[-1]
    for k in range(0,len(sing_w)-1):
        midw[k+1] = (sing_w[k] + sing_w[k+1])/2
    
    EE = np.zeros([Nc,len(sing_w)], dtype = 'complex128')
    
    for k in range(0,len(midw)):
        sk = 1j*midw[k]
        Anew = np.diag(A)
        Gnew = fitcalcABCDE(sk,Anew,B,C,D,E)
        G = np.real(Gnew)
        EE[:,k] = lin.eig(G)[0]
        if any(EE[:,k] < 0):
            viol[k] = 1
        else:
            viol[k] = 0
    
    intervals = np.zeros([2,1], dtype = 'double')
    for k in range(0,len(midw)):
        if viol[k] == 1:
            if k == 0:
                intervals[:,:] = np.transpose([0,sing_w[0]]).reshape(-1,1)
                #intervals = np.hstack((intervals,dummy))
            elif k == len(midw)-1:
                dummy = np.transpose([sing_w[k-1],1e16]).reshape(-1,1)
                intervals = np.hstack((intervals,dummy))
            else:
                dummy = np.transpose([sing_w[k-1],sing_w[k]]).reshape(-1,1)
                intervals = np.hstack((intervals,dummy))
    
    if len(intervals) == 0:
        wintervals = intervals
        return
    
    #Collapsing overlapping bands: DOUBLE CHECK THE BELOW
    tell = -1
    killindex = 0
    for k in range(1,0,len(intervals[0,:])):
        if intervals[1,k-1] == intervals[0,k]:
            tell = tell +1
            intervals[1,k-1] = intervals[1,k]
            intervals[:,k] = intervals[:,k-1]
            killindex[tell] = k - 1
    
    if killindex != 0:
        intervals[:,killindex] = []
    wintervals = intervals
    
    return wintervals
                