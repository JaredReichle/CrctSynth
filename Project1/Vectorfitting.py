"""
Vectorfitting.py
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

from pylab import zeros, vstack, concatenate, real, imag
from pylab import  outer, sort, ones, diag
import numpy as np
from numpy.linalg import eigvals, lstsq

#Complex conjugate function
def cc(z):
    return z.conjugate()

#Returns the transfer function
def model(s, poles, residues, d, h):
    return np.sum(r/(s-p) for p, r in zip(poles, residues)) + d + s*h

def vectfit_step(f, s, poles):
    """
    f = complex data to fit
    s = j*frequency
    poles = initial poles guess
        note: All complex poles must come in sequential complex conjugate pairs
    returns adjusted poles
    """
    N = len(poles)
    Ns = len(s)

    cindex = zeros(N)
    # cindex is:
    #   - 0 for real poles
    #   - 1 for the first of a complex-conjugate pair
    #   - 2 for the second of a cc pair
    for i, p in enumerate(poles):
        if p.imag != 0:
            if i == 0 or cindex[i-1] != 1:
                assert cc(poles[i]) == poles[i+1], ("Complex poles must come in conjugate pairs: %s, %s" % (poles[i], poles[i+1]))
                cindex[i] = 1
            else:
                cindex[i] = 2

    # First linear equation to solve. See Appendix A
    A = zeros((Ns, 2*N+2), dtype=np.complex64)
    for i, p in enumerate(poles):
        if cindex[i] == 0:
            A[:, i] = 1/(s - p)
        elif cindex[i] == 1:
            A[:, i] = 1/(s - p) + 1/(s - cc(p))
        elif cindex[i] == 2:
            A[:, i] = 1j/(s - p) - 1j/(s - cc(p))
        else:
            raise RuntimeError("cindex[%s] = %s" % (i, cindex[i]))

        A[:, N+2+i] = -A[:, i] * f[i]

    A[:, N] = 1
    A[:, N+1] = s

    # Solve Ax == b using pseudo-inverse
    b = f
    A = vstack((real(A), imag(A)))
    b = concatenate((real(b), imag(b)))
    x, residuals, rnk, s = lstsq(A, b, rcond=-1)

    residues = x[:N]
    d = x[N]
    h = x[N+1]

    # We only want the "tilde" part in (A.4)
    x = x[-N:]

    # Calculation of zeros: Appendix B
    A = diag(poles)
    b = ones(N)
    c = x
    for i, (ci, p) in enumerate(zip(cindex, poles)):
        if ci == 1:
            x, y = real(p), imag(p)
            A[i, i] = A[i+1, i+1] = x
            A[i, i+1] = -y
            A[i+1, i] = y
            b[i] = 2
            b[i+1] = 0
            #cv = c[i]
            #c[i,i+1] = real(cv), imag(cv)

    H = A - outer(b, c)
    H = real(H)
    new_poles = sort(eigvals(H))
    unstable = real(new_poles) > 0
    new_poles[unstable] -= 2*real(new_poles)[unstable]
    return new_poles

def calculate_residues(f, s, poles, rcond=-1):
    Ns = len(s)
    N = len(poles)

    cindex = zeros(N)
    for i, p in enumerate(poles):
        if p.imag != 0:
            if i == 0 or cindex[i-1] != 1:
                assert cc(poles[i]) == poles[i+1], ("Complex poles must come in conjugate pairs: %s, %s" % poles[i:i+1])
                cindex[i] = 1
            else:
                cindex[i] = 2

    # use the new poles to extract the residues
    A = zeros((Ns, N+2), dtype=np.complex128)
    for i, p in enumerate(poles):
        if cindex[i] == 0:
            A[:, i] = 1/(s - p)
        elif cindex[i] == 1:
            A[:, i] = 1/(s - p) + 1/(s - cc(p))
        elif cindex[i] == 2:
            A[:, i] = 1j/(s - p) - 1j/(s - cc(p))
        else:
            raise RuntimeError("cindex[%s] = %s" % (i, cindex[i]))

    A[:, N] = 1
    A[:, N+1] = s
    # Solve Ax == b using pseudo-inverse
    b = f
    A = vstack((real(A), imag(A)))
    b = concatenate((real(b), imag(b)))
    cA = np.linalg.cond(A)
    if cA > 1e13:
        print('Warning!: Ill Conditioned Matrix. Consider scaling the problem down')
        print('Cond(A)', cA)
    x, residuals, rnk, s = lstsq(A, b, rcond=rcond)

    # Recover complex values
    x = np.complex64(x)
    for i, ci in enumerate(cindex):
       if ci == 1:
           r1, r2 = x[i:i+2]
           x[i] = r1 - 1j*r2
           x[i+1] = r1 + 1j*r2

    residues = x[:N]
    d = x[N].real
    h = x[N+1].real
    return residues, d, h

def print_params(poles, residues, d, h):
    cfmt = "{0.real:g} + {0.imag:g}j"
    print("poles: " + ", ".join(cfmt.format(p) for p in poles))
    print("residues: " + ", ".join(cfmt.format(r) for r in residues))
    print("offset: {:g}".format(d))
    print("slope: {:g}".format(h))

def vectfit_auto(f, s, n_poles=10, n_iter=10, show=False,
                 inc_real=False, loss_ratio=1e-2, rcond=-1, track_poles=False):
    w = imag(s)
    pole_locs = np.linspace(w[0], w[-1], n_poles+2)[1:-1]
    lr = loss_ratio
    init_poles = poles = concatenate([[p*(-lr + 1j), p*(-lr - 1j)] for p in pole_locs])

    if inc_real:
        poles = concatenate((poles, [1]))

    poles_list = []
    for _ in range(n_iter):
        poles = vectfit_step(f, s, poles)
        poles_list.append(poles)

    residues, d, h = calculate_residues(f, s, poles, rcond=rcond)

    if track_poles:
        return poles, residues, d, h, np.array(poles_list)

    print_params(poles, residues, d, h)
    return poles, residues, d, h

def vectfit_auto_rescale(f, s, **kwargs):
    s_scale = abs(s[-1])
    f_scale = abs(f[-1])
    print('SCALED')
    poles_s, residues_s, d_s, h_s = vectfit_auto(f / f_scale, s / s_scale, **kwargs)
    poles = poles_s * s_scale
    residues = residues_s * f_scale * s_scale
    d = d_s * f_scale
    h = h_s * f_scale / s_scale
    print('UNSCALED')
    print_params(poles, residues, d, h)
    return poles, residues, d, h
