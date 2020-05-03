"""
Project3.py
Name: Jared Reichle
Date created: 4/12/2020
Description: This program will input a pole/residue arbitrary N-port network
that will then enforce passivity (eig(real(Y))>0) and output a passive (and stable)
N-port SPICE equivalent circuit that runs stable in time domain simulations.

Resources:
Code taken and modified from: https://www.sintef.no/projectweb/
vectorfitting/old-downloads/
Credit goes to the following papers:
    
    [1] A. Zadehgol, "A semi-analytic and cellular approach to rational system
    characterization through equivalent circuits", Wiley IJNM, 2015. [Online].
    https://doi.org/10.1002/jnm.2119

    [2]  V. Avula and A. Zadehgol, "A Novel Method for Equivalent Circuit
    Synthesis from Frequency Response of Multi-port Networks", EMC EUR, pp.
    79-84, 2016. [Online]. Available: ://WOS:000392194100012.
    
    [3] Gustavsen, B. (2008). "Fast Passivity Enforcement for Pole-Residue
    Models by Perturbation of Residue Matrix Eigenvalues." Power Delivery,
    IEEE Transactions on 23(4): 2278-2285.

    [4] Gustavsen, B. (2010). "Fast Passivity Enforcement for S-Parameter
    Models by Perturbation of Residue Matrix Eigenvalues." Advanced Packaging,
    IEEE Transactions on 33(1): 257-265.

    [5] Gustavsen, B. (2020), https://www.sintef.no/projectweb/vectorfitting/
"""
import numpy as np
from scipy import linalg as la
import matplotlib.pyplot as plt

SERA = np.genfromtxt('SERA.csv', delimiter=',', dtype = complex)
SERC = np.genfromtxt('SERC.csv', delimiter=',', dtype = complex)
SERC = np.reshape(SERC,(2,2,8), order = 'F')
SERD = np.genfromtxt('SERD.csv', delimiter=',', dtype = complex)
SERE = np.genfromtxt('SERE.csv', delimiter=',', dtype = complex)
#sPR = np.genfromtxt('sPR.csv', delimiter=',', dtype = complex)
sPR = 2*np.pi*1j*np.logspace(-1,7,301)

#Import bigY admittance array
#def BuildBigY():
Ns=501    #Number of frequency samples
Nc=2      #Size of Y (after reduction)
bigY = np.zeros([Nc,Nc,Ns], dtype = 'complex')
Y = np.zeros([4,4], dtype = 'complex')
s = 2*np.pi*1j*np.logspace(1,5,Ns)

#Component values: 
R1=1    
L1=1e-3 
C1=1e-6
R2=5    
L2=5e-3
R3=1    
C3=1e-6
L4=1e-3
R4=1e-2 
L5=20e-3
R6=10   
C6=10e-6
R7=1    
C7=2e-6

#Building Y, reduction:
for k in range(0,Ns):
    sk = s[k]
    y1 = 1/(R1+sk*L1+(1/(sk*C1)))
    y2=1/( R2+sk*L2 )
    y3=1/( R3+(1/(sk*C3)))
    y4=1/( R4+sk*L4)
    y5=1/(sk*L5)
    y6=1/( R6+(1/(sk*C6)))
    y7=1/( R7+(1/(sk*C7)))
    Y[0,0] = y1+y3
    Y[1,1] = y4
    Y[2,2] = y3 +y4 +y5 +y6
    Y[3,3] = y1 +y2 +y6 +y7
    Y[0,2] = -y3
    Y[0,3] = -y1
    Y[1,2] = -y4
    Y[2,0] = -y3
    Y[2,1] = -y4
    Y[2,3] = -y6
    Y[3,0] = -y1
    Y[3,2] = -y6
  
    #Eliminating nodes 3 and 4:
    Yred = Y[0:1,0:1] - Y[0:1,2:3]*Y[2:3,2:3]**(-1)*Y[2:3,0:1]
    bigY[:,:,k] = Yred
    
    #bigY[0,0,k] = bigY[0,0,k] - 1e4    #Modifying element [0,0] to check pass.


Nc = len(bigY[:,0,0])
NsPR = len(sPR)
N = len(SERA)
oldSERA = SERA
oldSERC = SERC
oldSERD = SERD
oldSERE = SERE

#=========================================
#Enforcing passivity by modifying SERD: 
#=========================================

Y = np.zeros([Nc,Nc], dtype = 'complex')
Dcorr = []
for k in range(0,NsPR):
    sk = sPR[k]
    for row in range(0,Nc):
        for col in range(0,Nc):
            Yfit = SERD[row,col] + sk*SERE[row,col]
            Yfit = Yfit + sum(np.squeeze(SERC[row,col,0:N])/(sk - SERA[0:N]))
            Y[row,col] = Yfit
    G = np.real(Y)
    [Edum,V] = la.eig(G)
    E = np.zeros([Nc,Nc], dtype = 'complex')
    for row in range(0,Nc):
        E[row,row] = Edum[row]
    Edum2 = E
    for row in range(0,Nc):
        if Edum2[row,row] > 0:
            Edum2[row,row] = 0
    Dcorr = -V*Edum2*la.inv(V)
    SERD = SERD + Dcorr

#=================================================
#Plotting eig(G), before and after passivity enf.:
#=================================================

EE = np.zeros([Nc,NsPR])
oldEE = np.zeros([Nc,NsPR])
oldY = np.zeros([Nc,Nc], dtype = 'complex')
for k in range(0,NsPR):
    sk = sPR[k]
    for row in range(0,Nc):
        for col in range(0,Nc):
            Yfit = SERD[row,col] + sk*SERE[row,col]
            Yfit = Yfit + sum(np.squeeze(SERC[row,col,0:N])/(sk-SERA[0:N]))
            Y[row,col] = Yfit
            Yfit = oldSERD[row,col] + sk*oldSERE[row,col]
            Yfit = Yfit + sum(np.squeeze(oldSERC[row,col,0:N])/(sk-oldSERA[0:N]))
            oldY[row,col] = Yfit
    G = np.real(Y)
    [Edum,V] = la.eig(G)
    E = np.zeros([Nc,Nc], dtype = 'complex')
    for row in range(0,Nc):
        E[row,row] = Edum[row]
    EE[:,k] = np.real(E[:,0])
    G = np.real(oldY)
    [Edum,V] = la.eig(G)
    E = np.zeros([Nc,Nc], dtype = 'complex')
    for row in range(0,Nc):
        E[row,row] = Edum[row]
    oldEE[:,k] = np.real(E[:,0])
    
plt.figure(figsize = (10,10), dpi=80, edgecolor ='k')
plt.semilogx(sPR/(2*np.pi*1j),np.real(np.transpose(oldEE)), '-r', label ='Before enforcement')
plt.semilogx(sPR/(2*np.pi*1j),np.real(np.transpose(EE)), '--b', label = 'After enforcement')
plt.legend()
plt.title('Eigenvalues of G vs. Frequency')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Eigenvalues of G')
plt.savefig('Eigenvalues.jpg', dpi=80)
plt.show()
    
#===============================
#Writting of the SPICE netlist:
#===============================  
    
"""
WRITING THE SPICE NETLIST
"""

TOLD = 1e-12
TOLE = 1e-12
TOLCA = 1e-12

cindex = np.zeros(N)
#Create spice netlist file
f = open("P3.sp","w+")

# KNOWN ERROR: PSPICE WILL NOT SUPPORT NESTED SUBCIRCUITS
for m in range(0,N):
    if np.imag(SERA[m]) != 0:
        if m == 1:
            cindex[m] = 1
        else:
            if cindex[m-1] == 0 or cindex[m-1] == 2:
                cindex[m] = 1
                cindex[m] = 2
            else:
                cindex[m] = 2

f.write('$VINTAGE,1/n')
f.write('C <BUS1><BUS2><BUS3><BUS4><  OHM       ><   milliH    ><   microF    >\n')
f.write('C \n')

for row in range(0,Nc):
    for col in range(0,Nc):
        strp = 'C  ({0},{1})\n'
        f.write(strp.format(row,col))
        if row == col: #Diagonal element
            dum = SERD[row,:]
            dum[col] = 0
            D = SERD[row,col]+sum(dum)
            if D>TOLD:
                R0 = 1/D
                f.write('  A____')
                strp = '{0}'
                f.write(strp.format(row))
                f.write('                   ')
                strp = '{0}\n'
                f.write(strp.format(R0))
            dum = SERE[row,:]
            dum[col]= 0 
            C0 = SERE[row,col] + sum(dum)
            if C0 > TOLE:
                f.write('  A____')
                strp = '{0}'
                f.write(strp.format(row))
                f.write('                                                  ')
                strp = '{0}'
                f.write(strp.format(1e6*C0))

            for m in range(0,N):
                M = m + (m<0)*2^32 #Convert to hex
                if cindex[m] == 0: #real pole
                    a1 = SERA[m]
                    dum = np.squeeze(SERC[row,:,m])
                    dum[col] = 0
                    c1 = SERC[row,col,m] + sum(dum)
                    L1 = 1/c1
                    R1 = -a1/c1

                    if (abs(c1/a1)) > TOLCA:
                        f.write('  A____')
                        strp = '{0}'
                        f.write(strp.format(row))
                        f.write('                   ')
                        strp = '{0}'
                        f.write(strp.format(R1))
                        strp = '{0}'
                        f.write(strp.format(1000*L1))
                elif cindex[m] == 1: #complex pole, 1st part
                    a1 = np.real(np.squeeze(SERA[m]))
                    a2 = np.imag(np.squeeze(SERA[m]))  
                    dum = np.squeeze(SERC[row,:,m])
                    dum[col]=0
                    dum = SERC[row,col,m] + sum(dum)
                    c1 = np.real(dum)
                    c2 = np.imag(dum)
                    L= 1/(2*c1)
                    dum = c1*a1+c2*a2
                    R=(-2*a1+2*dum*L)*L
                    C=(a1^2+a2^2+2*dum*R)*L
                    C=1/C
                    G = -2*dum*C*L
                    if (abs(c1/a1)) > TOLCA:
                        f.write('  A____')
                        strp = '{0}'
                        f.write(strp.format(row))
                        f.write('A')
                        strp = '{0}'
                        f.write(strp.format(M))
                        f.write('__')
                        strp = '{0}'
                        f.write(strp.format(row))
                        f.write('             ')
                        strp = '{0}'
                        f.write(strp.format(R))
                        strp = '{0}'
                        f.write(strp.format(1000*L))
                        f.write('  A')
                        strp = '{0}'
                        f.write(strp.format(M))
                        f.write('__')
                        strp = '{0}'
                        f.write(strp.format(row))
                        f.write('                   ')
                        strp = '{0}'
                        f.write(strp.format(1/G))
                        f.write('  A')
                        strp = '{0}'
                        f.write(strp.format(M))
                        f.write('__')
                        strp = '{0}'
                        f.write(strp.format(row))
                        f.write('                                                  ')
                        strp = '{0}'
                        f.write(strp.format(1e6*C))
    else: #row!=col (off-diagonal element)
        C0 = -SERE[row,col]
        if abs(SERD[row,col])>TOLD:
            R0=-1/SERD[row,col] 
            f.write('  A____')
            strp = '{0}' 
            f.write(strp.format(row))
            f.write('A____')
            strp = '{0}'
            f.write(strp.format(col))
            f.write('            ')
            strp = '{0}'
            f.write(strp.format(R0))       
        if abs(C0)>TOLE:
            f.write('  A____')
            strp = '{0}'
            f.write(strp.format(row))
            f.write('A____')
            strp = '{0}'
            f.write(strp.format(col))
            f.write('                                            ')
            strp = '{0}'
            f.write(strp.format(1e6*C0))
        for m in range(0,N):
            M = m+(m<0)*2^32 #Convert to hex
            if cindex[m] == 0: #real pole
                a1 = SERA[m]
                c1 = -SERC[row,col,m]
                L1=1/c1
                R1=-a1/c1
                if (abs(c1/a1))>TOLCA:
                    f.write('  A____')
                    strp = '{0}'
                    f.write(strp.format(row))
                    f.write('A____')
                    strp = '{0}'
                    f.write(strp.format(col))
                    f.write('            ')
                    strp = '{0}'
                    f.write(strp.format(R1))
                    strp = '{0}'
                    f.write(strp.format(1000*L1))         
            elif cindex[m] == 1: #complex pole, 1st part
                a1 = np.real(np.squeeze(SERA[m]))
                a2 = np.imag(np.squeeze(SERA[m]))  
                c1 = -np.real(SERC[row,col,m])
                c2 = -np.imag(SERC[row,col,m])
                L=2*c1
                L=1/L
                dum=c1*a1+c2*a2
                R=(-2*a1+2*dum*L)*L
                C=(a1^2+a2^2+2*dum*R)*L
                C=1/C
                G=-2*dum*C*L
                if (abs(c1/a1))>TOLCA:
                    f.write('  A____')
                    strp = '{0}'
                    f.write(strp.format(row))
                    f.write('A')
                    strp = '{0}'
                    f.write(strp.format(M))
                    f.write('_')
                    strp = '{0}{1}'
                    f.write(strp.format(row,col))
                    f.write('            ') 
                    strp = '{0}'
                    f.write(strp.format(R))
                    strp = '{0}'
                    f.write(strp.format(1000*L))
                    f.write('  A')
                    strp = '{0}'
                    f.write(strp.format(M))
                    f.write('_')
                    strp = '{0}{1}'
                    f.write(strp.format(row,col))
                    f.write('A____')
                    strp = '{0}'
                    f.write(strp.format(col))
                    f.write('            ')
                    strp = '{0}'
                    f.write(strp.format(1/G))
                    f.write('  A')
                    strp = '{0}'
                    f.write(strp.format(M))
                    f.write('_')
                    strp = '{0}{1}'
                    f.write(strp.format(row,col))
                    f.write('A____')
                    strp = '{0}'
                    f.write(strp.format(col))
                    f.write('                                            ')
                    strp = '{0}'
                    f.write(strp.format(1e6*C))
#Close the file
f.close()