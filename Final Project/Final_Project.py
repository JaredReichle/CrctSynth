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
    
    [4] A. Semlyen and B. Gustavsen, "A half-size singularity test matrix for
    fast and reliable passivity assessment of rational models", IEEE Trans.
    Power Del., vol. 24, no. 1, pp. 345-351, Jan. 2009.
    
a22F"""

from math import pi
import numpy as np
import scipy.io
from matplotlib import pyplot as plt

from VFDriver import VFDriver, DefVFOpts
from RPDriver import RPDriver, DefRPOpts
from NetListGen import NetListGen

#==============================
# IMPORT bigY AND s
#==============================

mat = scipy.io.loadmat('ex2_Y.mat')
bigY = mat['bigY']
s = mat['s']

#============================
#   VECTOR FITTING
#============================

poles = []
vfopts = DefVFOpts() #Imported from VFDriver
vfopts.N = 15

[SER, rmserr, bigYfit, opts2] = VFDriver(bigY, s, poles, vfopts)

#np.savez('data', R = SER.R, poles = SER.poles)
#=============================
#   PASSIVITY ENFORCEMENT
#=============================

rpopts = DefRPOpts() #Imported from RPDriver

[SER, bigYfit_passive, opts3] = RPDriver(SER, s, rpopts)

#=============================
#   NETLIST GENERATION
#=============================

filename = 'spice_net_list' #Do not add extension

#np.savez('data', R = SER.R, poles = SER.poles)

#Help with defining inputs to function
NetListGen(SER.R, SER.poles, filename)

#===================================================
#   COMPARING ORIGINAL MODEL WITH PERTURBED MODEL
#===================================================

plt.figure(8)
Nc = len(SER.D)
s = np.transpose(s)
for row in range(0,Nc):
    for col in range(row,Nc):
        dum0 = bigY[row,col,:]
        #dum1 = bigYfit[row,col,:]
        dum2 = bigYfit_passive[row,col,:]
        h1 = plt.semilogy(s/(2*pi*1j),abs(dum0),'b')
        h2 = plt.semilogy(s/(2*pi*1j),abs(dum2),'r--')
        h3 = plt.semilogy(s/(2*pi*1j),abs(dum2-dum0),'g-')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Admittance [S]')
plt.title('Original vs fitted model')
plt.legend((h1[0],h2[0],h3[0]),('Oringinal','Fitted','Difference'))
plt.savefig('images/Orivfitmodel.png', dpi = 80)

plt.figure(7)
Nc = len(SER.D)
#s = np.transpose(s)
for row in range(0,Nc):
    for col in range(row,Nc):
        #dum0 = bigY[row,col,:]
        dum1 = bigYfit[row,col,:]
        dum2 = bigYfit_passive[row,col,:]
        h1 = plt.semilogy(s/(2*pi*1j),abs(dum1),'b')
        h2 = plt.semilogy(s/(2*pi*1j),abs(dum2),'r--')
        h3 = plt.semilogy(s/(2*pi*1j),abs(dum2-dum1),'g-')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Admittance [S]')
plt.title('Fitted vs perturbated model')
plt.legend((h1[0],h2[0],h3[0]),('Fitted','Perturbated','Difference'))
plt.savefig('images/Fitvpermodel.png', dpi = 80)

