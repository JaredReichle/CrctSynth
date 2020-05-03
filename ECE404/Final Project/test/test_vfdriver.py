import pytest
from VFDriver import VFDriver,DefVFOpts
import numpy as np


def test_vfdriver():
    N=4
    opts = DefVFOpts()
    opts.screen = 1
    opts.N = 4
    bigH = np.zeros(shape=(N,N)) #(Nc,Nc,Ns) : matrix function to be fitted.
                                 #   Nc - demension of Hs
                                 #   Ns - number of frequency samples
    s = np.zeros(N)              # (1,Ns)  : vector of frequency points [rad/sec] 
    poles = np.zeros(N)          # (1,N)   :  vector of Initial poles. 

    R,a =  VFDriver(bigH, s, poles, opts)
        

