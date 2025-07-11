"""
ex4.py
Name: Jared Reichle
Date created: January 29th, 2020
Description: A four port system to be fitted to

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

test_s1 = 1j*np.linspace(1, 1e5, 800)
test_s2 = 1j*np.linspace(1, 1e5, 800)
test_s3 = 1j*np.linspace(1, 1e5, 800)
test_s4 = 1j*np.linspace(1, 1e5, 800)

test_poles1 = [
    -3100,
    -64000,
    -200+5000j, -200-5000j,
    -320+15000j, -320-15000j,
#    -4000+35000j, -4000-35000j,
#    -500+45000j, -500-45000j,
#    -6500+45000j, -6500-45000j,
#    -700+70000j, -700-70000j,
#    -8000+73000j, -8000-73000j,
#    -9000+90000j, -9000-90000j,
]
test_poles2 = [
    -2600,
    -76000,
    -100+5000j, -100-5000j,
    -120+15000j, -120-15000j,
#    -3000+35000j, -3000-35000j,
#    -200+45000j, -200-45000j,
#    -1500+45000j, -1500-45000j,
#    -500+70000j, -500-70000j,
#    -1000+73000j, -1000-73000j,
#    -2000+90000j, -2000-90000j,
]
test_poles3 = [
    -2100,
    -65000,
    -900+5000j, -900-5000j,
    -820+15000j, -820-15000j,
#    -7000+35000j, -7000-35000j,
#    -600+45000j, -600-45000j,
#    -5500+45000j, -5500-45000j,
#    -400+70000j, -400-70000j,
#    -3000+73000j, -3000-73000j,
#    -2000+90000j, -2000-90000j,
]
test_poles4 = [
    -7300,
    -13000,
    -100+5000j, -100-5000j,
    -320+15000j, -320-15000j,
#    -5000+35000j, -5000-35000j,
#    -700+45000j, -700-45000j,
#    -1500+45000j, -1500-45000j,
#    -500+70000j, -500-70000j,
#    -1000+73000j, -1000-73000j,
#    -2000+90000j, -2000-90000j,
]
test_residues1 = [
    -4300,
    -12000,
    -4+7000j, -4-7000j,
    -60+18000j, -60-18000j,
#    21000+45000j, 21000-45000j,
#    430+60000j, 430-60000j,
#    120+10000j, 120-10000j,
#    640000+80000j, 640000-80000j,
#    32000+45000j, 32000-45000j,
#    -2000+92000j, -2000-92000j
]
test_residues2 = [
    -5000,
    -64000,
    -5+7000j, -5-7000j,
    -60+18000j, -60-18000j,
#    32000+45000j, 32000-45000j,
#    10+60000j, 10-60000j,
#    650+10000j, 650-10000j,
#    20000+80000j, 20000-80000j,
#    3000+45000j, 3000-45000j,
#    -7000+92000j, -7000-92000j
]
test_residues3 = [
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
test_residues4 = [
    -98000,
    -54000,
    -5+7000j, -5-7000j,
    -540+18000j, -540-18000j,
#    32000+45000j, 32000-45000j,
#    670+60000j, 670-60000j,
#    340+10000j, 340-10000j,
#    10000+80000j, 10000-80000j,
#    5000+45000j, 5000-45000j,
#    -4000+92000j, -4000-92000j
]

test_d = .2
test_h = 2e-5

test_f1 = sum(c/(test_s1 - a) for c, a in zip(test_residues1, test_poles1))
test_f1 += test_d + test_h*test_s1
vectfit_auto(test_f1, test_s1)
poles1, residues1, d, h = vectfit_auto_rescale(test_f1, test_s1)
fitted1 = model(test_s1, poles1, residues1, d, h)

test_f2 = sum(c/(test_s2 - a) for c, a in zip(test_residues2, test_poles2))
test_f2 += test_d + test_h*test_s2
vectfit_auto(test_f2, test_s2)
poles2, residues2, d, h = vectfit_auto_rescale(test_f2, test_s2)
fitted2 = model(test_s2, poles2, residues2, d, h)

test_f3 = sum(c/(test_s3 - a) for c, a in zip(test_residues3, test_poles3))
test_f3 += test_d + test_h*test_s3
vectfit_auto(test_f3, test_s3)
poles3, residues3, d, h = vectfit_auto_rescale(test_f3, test_s3)
fitted3 = model(test_s3, poles3, residues3, d, h)

test_f4 = sum(c/(test_s4 - a) for c, a in zip(test_residues4, test_poles4))
test_f4 += test_d + test_h*test_s4
vectfit_auto(test_f4, test_s4)
poles4, residues4, d, h = vectfit_auto_rescale(test_f4, test_s4)
fitted4 = model(test_s4, poles4, residues4, d, h)

plt.figure()
plt.plot(test_s1.imag, test_f1.real)
plt.plot(test_s1.imag, test_f1.imag)
plt.plot(test_s1.imag, fitted1.real)
plt.plot(test_s1.imag, fitted1.imag)

plt.plot(test_s2.imag, test_f2.real)
plt.plot(test_s2.imag, test_f2.imag)
plt.plot(test_s2.imag, fitted2.real)
plt.plot(test_s2.imag, fitted2.imag)

plt.plot(test_s3.imag, test_f3.real)
plt.plot(test_s3.imag, test_f3.imag)
plt.plot(test_s3.imag, fitted3.real)
plt.plot(test_s3.imag, fitted3.imag)

plt.plot(test_s4.imag, test_f4.real)
plt.plot(test_s4.imag, test_f4.imag)
plt.plot(test_s4.imag, fitted4.real)
plt.plot(test_s4.imag, fitted4.imag)
plt.xlim([0, 30000])
plt.savefig('images/fourport.png')
plt.show()