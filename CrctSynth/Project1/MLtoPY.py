# ex1.m -> ex1.py            
#
# Fitting an artificially created frequency response (single element)
#
# -Creating a 3rd order frequency response f(s)
# -Fitting f(s) using vectfit3.m 
#   -Initial poles: 3 logarithmically spaced real poles -1 iteration
#
#This example script is part of the vector fitting package (VFIT3.zip) 
#Last revised: 08.08.2008. 
#Created by:   Bjorn Gustavsen and Jared Reichle

import numpy as np
import cmath
import Vectorfitting
from matplotlib import pyplot as plt

#Frequency samples:
Ns = 101
s = 2j*np.pi*np.logspace(0,4,Ns)
f = np.zeros(Ns)

print('Creating frequency response f(s)...') 
for k in range(Ns):
  sk = s[k]
  f[k] = 2/(sk+5) + (30+40j)/(sk-(-100+500j)) + (30-40j)/(sk-(-100-500j))+0.5


#Initial poles for Vector Fitting:
N=3 #order of approximation
poles = -2*np.pi*np.logspace(0,4,N) #Initial poles

weight = np.ones(Ns) #All frequency points are given equal weight

"""
opts.relax=1      #Use vector fitting with relaxed non-triviality constraint
opts.stable=1     #Enforce stable poles
opts.asymp=3      #Include both D, E in fitting    
opts.skip_pole=0  #Do NOT skip pole identification
opts.skip_res=0   #Do NOT skip identification of residues (C,D,E) 
opts.cmplx_ss=1   #Create complex state space model

opts.spy1=0       #No plotting for first stage of vector fitting
opts.spy2=1       #Create magnitude plot for fitting of f(s) 
opts.logx=1       #Use logarithmic abscissa axis
opts.logy=1       #Use logarithmic ordinate axis 
opts.errplot=1    #Include deviation in magnitude plot
opts.phaseplot=1  #Also produce plot of phase angle (in addition to magnitiude)
opts.legend=1     #Do include legends in plots
"""


print('vector fitting...')
#[SER,poles,rmserr,fit] = Vectorfitting.vectfit_auto(f,s,poles,weight) 
print('Done.')

poles, residues, d, h = Vectorfitting.vectfit_auto_rescale(f, s)
fitted = Vectorfitting.model(s, poles, residues, d, h)
plt.figure()
plt.xlim([0,2000])
plt.plot(s.imag, f.real)
plt.plot(s.imag, f.imag)
plt.plot(s.imag, fitted.real)
plt.plot(s.imag, fitted.imag)
plt.show()
