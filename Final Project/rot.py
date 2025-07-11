import numpy as np
from math import atan2, cos, sin, pi

def rot(S):
    Nc = len(S)
    SA = np.zeros([Nc,Nc], dtype = 'complex128')
    SB=SA
    scale1 = np.zeros(Nc, dtype = 'complex128')
    scale2 = scale1
    scale = scale1
    err1 = scale1
    err2 = scale1
    numerator = np.zeros(Nc)
    denominator = np.zeros(Nc)
    ang = np.zeros(Nc)
    
    for col in range(0,Nc):
        numerator[col] = 0.0
        denominator[col] = 0.0
        for j in range(0,Nc):
            numerator[col] = numerator[col]+np.imag(S[j,col])*np.real(S[j,col])
            denominator[col] = denominator[col]+(np.real(S[j,col]))**2 - (np.imag(S[j,col]))**2
        
        numerator[col] = -2.0*numerator[col]
        ang[col] = 0.5*atan2(numerator[col],denominator[col])
        
        scale1[col] = cos(ang[col]) + 1j*sin(ang[col])
        scale2[col] = cos(ang[col]+pi/2) + 1j*sin(ang[col]+pi/2)
        
        for j in range(0,Nc):
            SA[j,col] = S[j,col]*scale1[col]
            SB[j,col] = S[j,col]*scale2[col]
        
        #Square sum (error) of solution
        
        aaa = 0.0
        bbb = 0.0
        ccc = 0.0
        
        for j in range(0,Nc):
            aaa = aaa + (np.imag(SA[j,col]))**2
            bbb = bbb + np.real(SA[j,col])*np.imag(SA[j,col])
            ccc = ccc + (np.real(SA[j,col]))**2
            
        err1[col] = aaa*cos(ang[col])**2 + bbb*sin(2.0*ang[col]) + ccc*sin(ang[col])**2
        
        #Square sum (error) of solution #2:
        aaa = 0.0
        bbb = 0.0
        ccc = 0.0
        
        for j in range(0,Nc):
            aaa = aaa + (np.imag(SB[j,col]))**2
            bbb = bbb + np.real(SB[j,col])*np.imag(SB[j,col])
            ccc = ccc + (np.real(SB[j,col]))**2
            
        err2[col] = aaa*cos(ang[col])**2 + bbb*sin(2.0*ang[col]) + ccc*sin(ang[col])**2
        
        #Picking the solution (1,2) with the smallest square sum
        if(err1[col] < err2[col]):
            scale[col] = scale1[col]
        else:
            scale[col] = scale2[col]
            
        S[:,col] = S[:,col]*scale[col]
        
    OUT = S
    
    return OUT

