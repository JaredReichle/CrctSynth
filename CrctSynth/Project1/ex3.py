"""
ex3.py
Name: Jared Reichle
Date created: January 29th, 2020
Description: A single port system to be fitted to

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
import Vectorfitting
from matplotlib import pyplot as plt

test_s = 1j*np.linspace(1, 1e5, 800)

test_poles = [
    -4500,
    -41000,
    -100+5000j, -100-5000j,
    -120+15000j, -120-15000j,
#    -3000+35000j, -3000-35000j,
#    -200+45000j, -200-45000j,
#    -1500+45000j, -1500-45000j,
#    -500+70000j, -500-70000j,
#    -1000+73000j, -1000-73000j,
#    -2000+90000j, -2000-90000j,
]
test_residues = [
    -3000,
    -83000,
    -5+7000j, -5-7000j,
    -20+18000j, -20-18000j,
#    6000+45000j, 6000-45000j,
#    40+60000j, 40-60000j,
#    90+10000j, 90-10000j,
#    50000+80000j, 50000-80000j,
#    1000+45000j, 1000-45000j,
#    -5000+92000j, -5000-92000j
]
test_d = .2
test_h = 2e-5

test_f = sum(c/(test_s - a) for c, a in zip(test_residues, test_poles))
test_f += test_d + test_h*test_s
vectfit_auto(test_f, test_s)

poles, residues, d, h = vectfit_auto_rescale(test_f, test_s)
fitted = model(test_s, poles, residues, d, h)
plt.figure()
plt.xlim([0,30000])
plt.plot(test_s.imag, test_f.real)
plt.plot(test_s.imag, test_f.imag)
plt.plot(test_s.imag, fitted.real)
plt.plot(test_s.imag, fitted.imag)
plt.savefig('images/oneport.png')
plt.show()