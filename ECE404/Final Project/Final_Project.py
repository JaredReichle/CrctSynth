"""
Final_Project.py
Name: Jared Reichle
Date created: 4/30/2020
Description: 

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
from pylab import zeros, vstack, concatenate, real, imag
from pylab import  outer, sort, ones, diag
import numpy as np
from matplotlib import pyplot as plt

from VFDriver import VFDriver, DefVFOpts
from RPDriver import RPDriver, DefRPOpts

#TEST COMMENT #2

#==============================
#   POLE-RESIDUE FITTING
#==============================

class opts():
    N = 50
    poletype = 'linlogcmplx'
    weghtparam = 5
    Niter1 = 7
    Niter2 = 4
    asymp = 2
    logx = 0

poles = []
[SER, rmserr, bigYfit, opts2] = VFDriver(bigY, s, poles, opts)

#=============================
#   PASSIVITY ENFORCEMENT
#=============================

opts.parametertype = 'Y'
opts.plot.s_pass = np.transpose(2*pi*1j*np.linspace(0,2E5,1001))
opts.plot.ylim = [-2e-3, 2e-3]

[SER, bigYfit_passive, opts3] = RPDriver(SER, s, opts)

#=============================
#   WRITING TO NETLIST
#=============================

#Write code here from project 2


#===================================================
#   COMPARING ORIGINAL MODEL WITH PERTURBED MODEL
#===================================================
