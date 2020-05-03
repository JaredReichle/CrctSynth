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
N=3 #Order of approximation
poles = -2*np.pi*np.logspace(0,4,N) #Initial poles

weight = np.ones(Ns) #All frequency points are given equal weight


print('vector fitting...')

poles, residues, d, h = Vectorfitting.vectfit_auto_rescale(f, s)
fitted = Vectorfitting.model(s, poles, residues, d, h)

print('Done.')

plt.figure()
plt.xlim([0,2500])
plt.plot(s.imag, f.real)
plt.plot(s.imag, f.imag)
plt.plot(s.imag, fitted.real)
plt.plot(s.imag, fitted.imag)
plt.show()
