"""
Final_Project.py
Name: Jared Reichle
Date created: 4/30/2020
Description: See README.txt

Resources:
Credit for this program goes to Bjorn Gustaven.
Credit also goes to the following papers:
    
    [1] B. Gustavsen and A. Semlyen, "Rational approximation of frequency
    domain responses by Vector Fitting", IEEE Trans. Power Delivery,
    vol. 14, no. 3, pp. 1052-1061, July 1999.
    
    [2] A. Zadehgol, "A semi-analytic and cellular approach to rational
    system characterization through equivalent circuits", Wiley IJNM,
    2015. [Online]. https://doi.org/10.1002/jnm.2119

    [3] Gustavsen, B. (2008). "Fast Passivity Enforcement for Pole-Residue
    Models by Perturbation of Residue Matrix Eigenvalues." Power Delivery,
    IEEE Transactions on 23(4): 2278-2285.
"""

from math import pi
#from pylab import zeros, vstack, concatenate, real, imag
#from pylab import  outer, sort, ones, diag
import numpy as np
import scipy.io
#from matplotlib import pyplot as plt

from VFDriver import VFDriver, DefVFOpts
from RPDriver import RPDriver, DefRPOpts
from NetListGen import NetListGen

#==============================
# IMPORT bigY AND s
#==============================

mat = scipy.io.loadmat('ex2_Y.mat')
bigY = mat['bigY']
s = mat['s']

#class opts():
#    N = 10
#    nu = 1e-3
#    poletype = 'linlogcmplx'
#    weightparam = 1
#    Niter1 = 7
#    Niter2 = 4
#    asymp = 2
#    logx = 0
#    remove_HFpoles = 0
#    factor_HF = 1.1
#    passive_DE = 0
#    passive_DE_TOLD = 1e-6
#    passive_DE_TOLE = 1e-16
#    screen = 1

#============================
#   VECTOR FITTING
#============================

poles = []
vfopts = DefVFOpts()

[SER, rmserr, bigYfit, opts2] = VFDriver(bigY, s, poles, vfopts)

np.savez('data', R = SER.R, poles = SER.poles)
#=============================
#   PASSIVITY ENFORCEMENT
#=============================

filename = 'spice_net_list' #Do not add extension

#Help with defining inputs to function
NetListGen(SER.R, SER.poles, filename)


rpopts = DefRPOpts()

#opts.parametertype = 'Y'
#opts.plot.s_pass = np.transpose(2*pi*1j*np.linspace(0,2E5,1001))
#opts.plot.ylim = [-2e-3, 2e-3]

[SER, bigYfit_passive, opts3] = RPDriver(SER, s, rpopts)
#=============================
#   NETLIST GENERATION
#=============================

#===================================================
#   COMPARING ORIGINAL MODEL WITH PERTURBED MODEL
#===================================================

plt.figure(1)
Nc = len(SER.D)
for row in range(0,Nc):
    for col in range(row,Nc):
        dum1 = np.squeeze(bigYfit[row,col,:])
        dum2 = np.squeeze(bigYfit_passive[row,col,:])
        h1 = plt.semilogy(s/(2*pi*1j),abs(dum1),'b')
        h2 = plt.semilogy(s/(2*pi*1j),abs(dum2),'r--')
        h3 = plt.semilogy(s/(2*pi*1j),abs(dum2-dum1),'g-')

plt.xlabel('Frequency [Hz]')
plt.ylabel('Admittance [S]')
plt.legend()
