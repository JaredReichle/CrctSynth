import numpy as np

def NetListGen(residues, poles, filename):
    
    Np = len(poles)
    Nr = len(residues[0,0,:])
    
    PORT = len(residues[0])
    
    if Np != Nr:
        print('ERROR in NetListGen.py: poles is not the same length as residues')
        return
    else:
        NN = Np
    
    #Initializing P/R
    pi = np.zeros(NN)
    pr = np.zeros(NN)
    ci = np.zeros([PORT,NN])
    cr = np.zeros([PORT,NN])
    
    
    for n in range(NN): #Poles
        pi[n] = np.imag(poles[n])
        pr[n] = np.real(poles[n])
    
    for m in range(0,PORT): #Residues
        for l in range(NN):
            ci[m,l] = np.imag(residues[0,m,l])
            cr[m,l] = np.real(residues[0,m,l])    
    
    #Initializing RLRC vars
    Racplx = np.zeros([PORT,NN])
    Rbcplx = np.zeros([PORT,NN])
    Ccplx = np.zeros([PORT,NN])
    Lcplx = np.zeros([PORT,NN])
    
    #Initializing RL vars
    Rreal = np.zeros([PORT,NN])
    Lreal = np.zeros([PORT,NN])
    
    """
    POLE AND RESIDUE RLRC AND RL RELATIONS
    """
    
    ### FOR REAL P/R ###
    #R = -pr/cr
    #L = 1/cr
    
    ### FOR COMPLEX CONJUGATE P/R ###
    #Ra = (ci*pi - cr*pr)/(2*cr**2)
    #Rb = -(pi**2*(ci**2+cr**2))/(2*cr**2*(ci*pi+cr*pr))
    #L = 1/(2*cr)
    #C = (2*cr**3)/(pi**2*(ci**2+cr**2))
    
    #Sanity checks - for better optimization, implement Table III from [2]
    for p in range(0,PORT):
        for n in range(0,NN):
            #Check stability, causality, and passivity
            if (pr[n] > 0):
                print('ERROR in NetListGen.py: Accessed pole is not stable or causal. Returning...')
                print('Error found at pole number: ')
                print(n)
                return
            if (pr[n] < 0 and cr[p,n] < 0):
                print('ERROR in NetListGen.py: Accessed pole/residue pair is not passive. Returning...')
                print('Error found at pole number: ')
                print(n)
                print('and port number: ')
                print(p)
                return
            
    #Passed above checks, Assign values:
    for p in range(0,PORT):
        for n in range(0,NN):
            if(pi[n] == 0 and ci[p,n] == 0 and cr[p,n]!=0 and pr[n]!=0): #RL BRANCH
                Rreal[p,n] = -pr[n]/cr[p,n]
                Lreal[p,n] = 1/cr[p,n]
            else: #RLRC BRANCH
                if(cr[p,n]==0 and ci[p,n]==0 and pr[n]==0 and pi[n]==0): #A non-sense branch
                    print('ERROR in NetListGen.py: branch creation failed')
                    return
                else: #abs added for simulation reasons
                    Racplx[p,n] = abs((ci[p,n]*pi[n] - cr[p,n]*pr[n])/(2*(cr[p,n]**2)))
                    Rbcplx[p,n] = abs(((pi[n]**2)*((ci[p,n]**2)+(cr[p,n]**2)))/(2*(cr[p,n]**2)*(ci[p,n]*pi[n]+cr[p,n]*pr[n])))
                    Lcplx[p,n] = abs(1/(2*cr[p,n]))
                    Ccplx[p,n] = abs((2*(cr[p,n]**3))/(pi[n]**2*(ci[p,n]**2+cr[p,n]**2)))
    
    """
    WRITING THE SPICE NETLIST
    """
    
    #Create spice netlist file
    fileext = '.sp'
    fileopen = str(filename)
    file = fileopen+fileext
    f = open(file,"w+")
    
    #Subcircuit file
    # KNOWN ERROR: PSPICE WILL NOT SUPPORT NESTED SUBCIRCUITS
    # IF DOING IN PSPICE, COMMENT BELOW TWO LINES AND NESTED .ENDS
    #sub = (".SUBCKT PORT NODE_TOP NODE_BOT\n")
    #f.write(sub.format(0))
    
    # Writing all blocks (ports) for system
    for u in range(1,PORT+1):
        # Writing all branches within block
        Subckt = ".SUBCKT PORT_{0} NODE_TOP NODE_BOT\n"
        f.write(Subckt.format(u))
        for i in range(0,NN):
            if Rreal[u-1,i] != 0 and Lreal[u-1,i] != 0: #Write the RL branch
                Rrealstr = "R{0} NODE_TOP NODE_MID {1}\n"
                Lrealstr = "L{0} NODE_MID NODE_BOT {1}\n"
                f.write(Rrealstr.format(i,Rreal[i]))
                f.write(Lrealstr.format(i,Lreal[i]))
            else:                                       #Write the RLRC branch
                Rastr = "RA{0} NODE_TOP NODE_{1}2 {2}\n"
                Rbstr = "RB{0} NODE_{1}3 NODE_BOT {2}\n"
                Lstr = "L{0} NODE_{1}2 NODE_{1}3 {2}\n"
                Cstr = "C{0} NODE_{1}3 NODE_BOT {2}\n"
                f.write(Rastr.format(i, i, Racplx[u-1,i]))
                f.write(Rbstr.format(i, i, Rbcplx[u-1,i]))
                f.write(Lstr.format(i, i, Lcplx[u-1,i]))
                f.write(Cstr.format(i, i, Ccplx[u-1,i]))
        #f.write(".PLOT AC VM(NODE_TOP)\n")
        ends = ".ENDS PORT_{0}\n\n"
        f.write(ends.format(u))
    #END MAIN SUBCIRCUIT - COMMENT IF USING PSPICE
    #f.write(".ENDS\n")
    #Termination
    f.write(".END")
    #Close the file
    f.close()
    return
    
def main(argv):
    ofilename = argv[1]
    ifilename = argv[2]
    data = np.load(ifilename)
    R = data['R']
    poles = data['poles']
    NetListGen(R, poles, ofilename)
    
if __name__ == '__main__':
    import sys
    main(sys.argv)