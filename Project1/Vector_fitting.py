"""
Vector_fitting.py
Name: Jared Reichle
Date created: January 23rd, 2020
Description: Approximate f(s) using the vector fitting algorithm

Resources:
Credit for this program goes to Bjorn Gustaven and Phil Reinhold.
Code taken and modified from: https://github.com/PhilReinhold/vectfit_python/blob/
master/vectfit.py
Credit also goes to the following papers:
    
    [1] B. Gustavsen and A. Semlyen, "Rational approximation of frequency
    domain responses by Vector Fitting", IEEE Trans. Power Delivery,
    vol. 14, no. 3, pp. 1052-1061, July 1999.
    
    [2] B. Gustavsen, "Improving the pole relocating properties of vector 
    fitting", IEEE Trans.Power Delivery, vol 21, no 3, pp 1587-1592, July 2006
    
    [3] D. Deschrijver, M. Mrozowski, T. Dhaene and D. de Zutter: "Macromodeling 
    of Multiport Systems Using a Fast Implementation of the Vector Fitting Method", 
    IEEE Microwave and Wireless Components Letters, vol. 18, no 6, pp 383-385, 
    June 2008
    
    [4] T. D. Nguyen, “Stability and Passivity Analysis, Application of Vector
    Fitting in RLC Equivalent Circuit Synthesis,” M.S. thesis, ECE, UIdaho,
    Moscow, ID, 2016.
"""

import numpy as np

""" Objects """
Ns = 101 #Number of frequency samples
f = np.mgrid[1:1e5:Ns, 1:1e5:Ns].reshape(2,-1).T
N = 10 # Order of approx
w = f*(2*np.pi)
s = 1j*w


bet = np.linspace(w[1], w[Ns-1], N/2)
m = len(bet)
poles = np.zeros((m,2))
for n in range(m):
    alf = -bet[n]*1e-2
    poles[n][0] = (alf-1j*bet[n])
    poles[n][1] = (alf+1j*bet[n])

weight = np.ones([Ns,Ns])

class opts():
    relax = 1       #Use vector fitting with relaxed non-triviality constraint
    stable=1        #Enforce stable poles
    asymp=2         #Include only D in fitting (not E)   
    skip_pole=0     #Do NOT skip pole identification
    skip_res=0      #Do NOT skip identification of residues (C,D,E) 
    cmplx_ss=1      #Create complex state space model
    spy1=0          #No plotting for first stage of vector fitting
    spy2=1          #Create magnitude plot for fitting of f(s) 
    logx=1          #Use logarithmic abscissa axis
    logy=1          #Use logarithmic ordinate axis 
    errplot=1       #Include deviation in magnitude plot
    phaseplot=0     #exclude plot of phase angle (in addition to magnitiude)
    legend=1        #Do include legends in plots

TOLlow = 1e-18
TOLhigh = 1e18

a, b = poles.shape
if s[0, 0] == 0 and a == 1:
    if poles[0] == 0 and poles[1] != 0:
        poles[0] = -1
    elif poles[1] == 0 and poles[0] != 0:
        poles[1] = -1
    elif poles[1] == 0 and poles[1] == 0:
        poles[0] = -1 + 10j
        poles[1] = -1 - 10j

rmserr = (np.array([])).reshape(0, 0)
a, b = s.shape
if a < b:
    s = s.T

LAMBD = np.diagflat(poles)
Ns = len(s)
N = len(LAMBD)
Nc = len(f[:, 0])
B = np.ones([N, 1])
SERA = poles.copy()
SERC = np.zeros([Nc, N])
SERD = n
p.zeros([Nc, 1])
SERE = np.zeros([Nc, 1])
roetter = poles.copy()
fit = np.zeros([Nc, Ns])

weight = weight.T
if len(weight[0, :]) == 1:
    common_weight = 1
elif len(weight[1, :]) == Nc:
    common_weight = 0

if opts.asymp == 1:
    offs = 0
elif opts.asymp == 2:
    offs = 1
else:
    offs = 2

# Pole Identification
Escale = np.zeros([1, Nc+1])
cindex = np.zeros([1, N])
    # determine complex poles
for m in range(N):
    if np.imag(LAMBD[m, m]) != 0:
        if m == 0:
            cindex[m] = 1
        else:
            if cindex[m-1] == 0 or cindex[m-1] == 2:
                cindex[m] = 1
                cindex[m+1] = 2
            else:
                cindex[m] = 2

    # build system matrix
Dk = np.zeros([Ns, N], dtype='complex')
for m in range(N):
    if cindex[0, m] == 0:
        Dk[:, m] = (1/(s-LAMBD[m, m])).reshape(len(s),)
    elif cindex[0, m] == 1:
        Dk[:, m] = 1/(s-LAMBD[m, m]) + 1/(s-np.conj(LAMBD[m, m]))
        Dk[:, m+1] = 1j/(s-LAMBD[m, m]) - 1j/(s-np.conj(LAMBD[m, m]))

if opts.asymp == 1 or opts.asymp == 2:
    Dk = np.hstack((Dk, np.ones([len(Dk), 1])))
elif opts.asymp == 3:
    Dk = np.hstack((Dk, np.ones([len(Dk), 1])))
    Dk = np.hstack((Dk, s))

scale = 0
for m in range(Nc):
    if len(weight[0, :]) == 1:
        scale = scale+(np.linalg.norm(weight[:, 0]*f[m, :].T))**2
    else:
        scale = scale+(np.linalg.norm(weight[:, m]*f[m, :].T))**2
scale = np.sqrt(scale)/Ns

# # relaxed part

AA = np.zeros([Nc*(N+1), N+1])
bb = np.zeros([Nc*(N+1), 1])
Escale = np.zeros([1, len(AA[0, :])])
for n in range(Nc):
    A = np.zeros([Ns, (N+offs)+N+1], dtype='complex')
    if common_weight == 1:
        weig = weight
    else:
        weig = weight[:, n]
    for m in range(N+offs):
        A[0:Ns, m] = weig[m]*Dk[0:Ns, m].T
    inda = N+offs
    for m in range(N+1):
        A[0:Ns, inda+m] = -weig[m]*Dk[0:Ns, m]*f[n, 0:Ns].T
    A = np.vstack((np.real(A), np.imag(A)))
    offset = N+offs
    if n == Nc-1:
        ap = np.zeros(len(A[0, :]))
        for mm in range(N+1):
            ap[offset+mm] = np.real(scale*np.sum(Dk[:, mm]))
        A = np.vstack((A, ap))
    Q, R = np.linalg.qr(A)
    ind1 = N+offs
    ind2 = N+offs+N+1
    R22 = R[ind1:ind2, ind1:ind2]
    AA = R22.copy()
    if n == Nc-1:
        for k in range(len(Q[0,:])-(N+offs)):
            bb[k,0] = Q[-1, k+N+offs]*Ns*scale

for col in range(len(AA[0, :])):
    Escale[0, col] = 1/np.linalg.norm(AA[:, col])
    AA[:, col] = Escale[0, col]*AA[:, col]
x = np.linalg.solve(AA, bb)
x = x*np.conj(Escale.T)

C = x[0:-1]
D = x[-1]
for m in range(N):
    if cindex[0, m] == 1:
        r1 = C[m]
        r2 = C[m+1]
        C[m] = r1 + 1j*r2
        C[m+1] = r1 - 1j*r2

m = 0
for n in range(N):
    if m < N:
        if abs(LAMBD[m, m]) > abs(np.real(LAMBD[m, m])):
            LAMBD[m+1, m] = -np.imag(LAMBD[m, m])
            LAMBD[m, m+1] = np.imag(LAMBD[m, m])
            LAMBD[m, m] = np.real(LAMBD[m, m])
            LAMBD[m+1, m+1] = LAMBD[m, m]
            B[m, 0] = 2
            B[m+1, 0] = 0
            koko = C[m]
            C[m] = np.real(koko)
            C[m+1] = np.imag(koko)
            m += 1
    m += 1
B = np.matrix(B)
C = np.matrix(C)
D = np.matrix(D)
ZER = LAMBD - B*C.T/D
B = np.array(B)
C = np.array(C)
D = np.array(D)
ZER = np.array(ZER)
roetter, hold = np.linalg.eig(ZER)
roetter = roetter.reshape(1, 3)
unstables = np.real(roetter)>0
if opts.stable == 1:
    roetter[unstables] = roetter[unstables] - 2*np.real(roetter[unstables])
roetter = np.sort(roetter)
N = len(roetter[0,:])
for n in range(N):
    for m in range(n+1, N):
        if np.imag(roetter[0, m]) == 0 and np.imag(roetter[0, n]) != 0:
            trans = roetter[0, n]
            roetter[0, n] = roetter[0, m]
            roetter[0, m] = trans
N1 = 0
for m in range(N):
    if np.imag(roetter[0, m]) == 0:
        N1 = m

if N1 < N:
    roetter[N1:N+1] = np.sort(roetter[N1:N+1])

SERA = np.flip(roetter.T)

# Residue Identification

LAMBD = np.flip(roetter)
cindex = np.zeros([1, N])
for m in range(N):
    if np.imag(LAMBD[0, m]) != 0:
        if m == 0:
            cindex[m] = 1
        else:
            if cindex[0, m-1] == 0 or cindex[0, m-1] == 2:
                cindex[0, m] = 1
                cindex[0, m+1] = 2
            else:
                cindex[0, m] = 2

if opts.asymp == 1:
    A = np.zeros([Ns, N])
    BB = np.zeros([Ns, Nc])
elif opts.asymp == 2:
    A = np.zeros([Ns, N+1])
    BB = np.zeros([Ns, Nc])
else:
    A = np.zeros([Ns, N+2])
    BB = np.zeros([Ns, Nc])

Dk = np.zeros([Ns, N], dtype='complex')
for m in range(N):
    if cindex[0, m] == 0:
        Dk[:, m] = (1/(s-LAMBD[0, m])).reshape(Ns,)
    elif cindex[0, m] == 1:
        Dk[:, m] = (1/(s-LAMBD[0, m]) + 1/(s-np.conj(LAMBD[0, m]))).reshape(Ns,)
        Dk[:, m+1] = (1j/(s-LAMBD[0, m]) - 1j/(s-np.conj(LAMBD[0, m]))).reshape(Ns,)

if common_weight == 1:
    Dk = weight*Dk

    A = A.astype('complex')
    if opts.asymp == 1:
        A[0:Ns, 0:N] = Dk
    elif opts.asymp == 2:
        A[0:Ns, 0:N] = Dk
        A[0:Ns, N+1] = weight
    else:
        A[0:Ns, 0:N] = Dk
        A[0:Ns, N] = weight.reshape(Ns,)
        A[0:Ns, N+1] = (weight*s).reshape(Ns,)

    BB = BB.astype('complex')
    for m in range(Nc):
        BB[0:Ns, m] = (weight*(f.reshape(Ns,1))).reshape(Ns,)

    A = np.vstack((np.real(A), np.imag(A)))
    BB = np.vstack((np.real(BB), np.imag(BB)))

    Escale = np.zeros([1, len(A[0,:])])
    for col in range(len(A[0,:])):
        Escale[0, col] = np.linalg.norm(A[:, col])
        A[:, col] = A[:, col]/Escale[0, col]

    X = np.linalg.lstsq(A, BB, rcond=None)
    X = (X[0]/Escale.T).T

    C = X[:, 0:N]
    if opts.asymp == 2:
        SERD = X[:, N]
    elif opts.asymp == 3:
        SERE = X[:, N+1]
        SERD = X[:, N]

else:
    print('use common weight')

C = C.astype('complex')
for m in range(N):
    if cindex[0, m] == 1:
        for n in range(Nc):
            r1 = C[n, m]
            r2 = C[n, m+1]
            C[n, m] = r1 + 1j*r2
            C[n, m+1] = r1 - 1j*r2

B = np.ones([N, 1])
SERA = LAMBD
SERB = B
SERC = C

Dk = np.zeros([Ns, N], dtype='complex')
for m in range(N):
    Dk[:, m] = (1/(s-SERA[0, m])).reshape(Ns,)

SERC = SERC.reshape(3, 1)
fit = fit.astype('complex')
for n in range(Nc):
    fit[n, :] = ((Dk.dot(SERC.T[n, :])).reshape(Ns, 1)).T
    if opts.asymp == 2:
        fit[n,:] = fit[n,:] + SERD[n]
    elif opts.asymp == 3:
        fit[n, :] = fit[n, :] + SERD[n] + (s.T)*SERE[n]

fit = fit.T
f = f.T
diff = fit-f
rmserr = np.sqrt(np.sum(np.sum(abs(diff**2))))/np.sqrt(Nc*Ns)
# fit = fit.T

A = SERA
poles = A.copy()
if opts.skip_res != 1:
    B = SERB
    C = SERC
    D = SERD
    E = SERE
else:
    B = np.ones([N, 1])
    C = np.zeros([Nc, N])
    D = np.zeros([Nc, Nc])
    E = np.zeros([Nc, Nc])
    rmserr = 0
