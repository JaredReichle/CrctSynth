import pytest
#from RPDriver import RPDriver,DefRPOpts
from RPDriver import RPDriver
import numpy as np

class DefRPOpts:
    N = 0
    screen = 1


def test_rpdriver():
    N=4
    opts = DefRPOpts()
    opts.screen = 1
    opts.N = 4
    class SER:
        poles = None
        R = None
        D = None
        E = None 

    mySER = SER()
        
    SER.poles = np.zeros(N) # (1,N)     : poles 
    SER.R = np.zeros(N,N)   # (Nc,Nc,N) : residue 
    SER.D = np.zeros(N,N)     # (Nc,Nc)   : constant term  
    SER.E = np.zeros(N,N)     # (Nc,Nc)   : capacitive term 

    s = np.zeros(N)         # (1,Ns)    : frequency samples (jw [rad/sec]) 
    poles = np.zeros(N)

    newSER,opts2 =  RPDriver( SER, s, opts)
        

