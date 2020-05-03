# ex2.m -> ex2.py        
#
# -Creating an 18th order frequency response f(s) of 2 elements.
# -Fitting f(s) using vectfit3.m 
#   -Initial poles: 9 linearly spaced complex pairs (N=18)
#   -3 iterations
#
# This example script is part of the vector fitting package (VFIT3.zip) 
# Last revised: 08.08.2008. 
# Created by:   Bjorn Gustavsen.
#

import numpy as np
import Vectorfitting
from matplotlib import pyplot as plt

D=0.2
E=2e-5
p=[-4500,  -41000, (-100+5e3j), (-100-5e3j), (-120+15e3j), (-120-15e3j), (-3e3+35e3j),
(-3e3-35e3j), (-200+45e3j), (-200-45e3j), (-1500+45e3j), (-1500-45e3j), (-5e2+70e3j),
(-5e2-70e3j), (-1e3+73e3j), (-1e3-73e3j),  (-2e3+90e3j), (-2e3-90e3j)]
r=[-3000, -83000, (-5+7e3j), (-5-7e3j), (-20+18e3j), (-20-18e3j), (6e3+45e3j),
(6e3-45e3j), (40 +60e3j), (40-60e3j), (90 +10e3j), (90-10e3j), (5e4+80e3j),
(5e4-80e3j), (1e3+45e3j), (1e3-45e3j), (-5e3+92e3j), (-5e3-92e3j)]
#p=[p, (-200+45e3j), (-200-45e3j), (-1500+45e3j), (-1500-45e3j)]
#r=[r, (40 +60e3j), (40-60e3j), (90 +10e3j), (90-10e3j)]
#p=[p, (-5e2+70e3j), (-5e2-70e3j), (-1e3+73e3j), (-1e3-73e3j)  (-2e3+90e3j) (-2e3-90e3j)]
#r=[r, (5e4+80e3j), (5e4-80e3j), (1e3+45e3j), (1e3-45e3j), (-5e3+92e3j), (-5e3-92e3j)]      

w = 2*np.pi*np.linspace(1,1e5,100)
Ns = len(w)
s = w*1j

p = [i * (2*np.pi) for i in p]
r = [j * (2*np.pi) for j in r]

p1 = p[1:10]
r1 = r[1:10]
N1 = len(p1)
p2 = p[9:18]
r2 = r[9:18]
N2 = len(p2)

f = np.zeros(Ns)
for k in range(Ns):
    for n in range(N1):
        f[k]=f[k]+r1[n]/(s[k]-p1[n])  
    f[k]=f[k]+s[k]*E
f[:] = f[:] + D

for k in range(Ns):
    for n in range(N2):
        f[k]=f[k]+r2[n]/(s[k]-p2[n])
    f[k]=f[k]+s[k]*3*E
f[:] = f[:] + 2 * D

#=====================================
# Rational function approximation of f(s):
#=====================================


N=18 #Order of approximation 

#Complex starting poles :
"""
bet = np.linspace(w[1],w[Ns-1],N/2)
poles = np
l = len(bet)
for n in range(l):
    alf = -bet[n]*1e-2
    poles = [poles (alf-bet(n)*1j) (alf+bet(n)*1j) ]
"""
# Real starting poles :
#poles=-linspace(w(1),w(Ns),N); 
 
#Parameters for Vector Fitting : 

weight = np.ones(Ns)

print('vector fitting...')

poles, residues, d, h = Vectorfitting.vectfit_auto_rescale(f, s)
fitted = Vectorfitting.model(s, poles, residues, d, h)

print('Done.')

plt.figure()
#plt.xlim([0,10000])
plt.plot(s.imag, f.real)
plt.plot(s.imag, f.imag)
plt.plot(s.imag, fitted.real)
plt.plot(s.imag, fitted.imag)
plt.show()