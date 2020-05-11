import numpy as np

def NetListGen(residues, poles, filename):
    
    #Initializing P/R
    
    pi = np.zeros(NN)
    pr = np.zeros(NN)
    ci = np.zeros([PORT-1,NN])
    cr = np.zeros([PORT-1,NN])
    
    
    for n in range(NN): #Poles
        pi[n] = np.imag(df[0,n])
        pr[n] = np.real(df[0,n])
    
    for m in range(1,PORT): #Residues
        for l in range(NN):
            ci[m-1,l] = np.imag(df[m,l])
            cr[m-1,l] = np.real(df[m,l])    
    
    #Initializing RLRC vars
    Ra = np.zeros([PORT-1,NN])
    Rb = np.zeros([PORT-1,NN])
    C = np.zeros([PORT-1,NN])
    L = np.zeros([PORT-1,NN])
    
    
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
    
    #Below is the calculations and decision process
    
    PORT = PORT - 1
    
    for p in range(0,PORT):
        for n in range(NN):
            if(pi[n] == 0 and ci[p,n] == 0 and cr[p,n]!=0 and pr[n]!=0): #RL BRANCH
                Ra[p,n] = -pr[n]/cr[p,n]
                Rb[p,n] = 0
                C[p,n] = 0
                L[p,n] = 1/cr[p,n]
            else: #RLRC BRANCH
                if(cr[p,n]==0 and ci[p,n]==0 and pr[n]==0 and pi[n]==0): #A non-sense branch
                    Ra[p,n] = 0
                    Rb[p,n] = 0
                    C[p,n] = 0
                    L[p,n] = 0
                else:
                    Ra[p,n] = (ci[p,n]*pi[n] - cr[p,n]*pr[n])/(2*(cr[p,n]**2))
                    Rb[p,n] = -((pi[n]**2)*((ci[p,n]**2)+(cr[p,n]**2)))/(2*(cr[p,n]**2)*(ci[p,n]*pi[n]+cr[p,n]*pr[n]))
                    if cr[p,n] == 0:
                        L[p,n] = 0
                    else:    
                        L[p,n] = 1/(2*cr[p,n])
                    if pi[n] == 0:
                        C[p,n] = 0
                    else:
                        C[p,n] = (2*(cr[p,n]**3))/(pi[n]**2*(ci[p,n]**2+cr[p,n]**2))  
    
    
    """
    WRITING THE SPICE NETLIST
    """
    #Create spice netlist file
    f = open("spice_net_list.sp","w+")
    
    #Subcircuit file
    # KNOWN ERROR: PSPICE WILL NOT SUPPORT NESTED SUBCIRCUITS
    # IF DOING IN PSPICE, COMMENT BELOW TWO LINES AND NESTED .ENDS
    #sub = (".SUBCKT PORT NODE_TOP NODE_BOT\n")
    #f.write(sub.format(0))
    
    # Writing all blocks (ports) for system
    for u in range(PORT):
        # Writing all branches within block
        Subckt = ".SUBCKT PORT_{0} NODE_TOP NODE_BOT\n"
        f.write(Subckt.format(u))
        for i in range(NN):
            Rastr = "RA{0} NODE_TOP NODE_{1}2 {2}\n"
            Rbstr = "RB{0} NODE_{1}3 NODE_BOT {2}\n"
            Lstr = "L{0} NODE_{1}2 NODE_{1}3 {2}\n"
            Cstr = "C{0} NODE_{1}3 NODE_BOT {2}\n"
            f.write(Rastr.format(i, i, Ra[u,i]))
            f.write(Rbstr.format(i, i, Rb[u,i]))
            f.write(Lstr.format(i, i, L[u,i]))
            f.write(Cstr.format(i, i, C[u,i]))
        f.write(".PLOT AC VM(NODE_TOP)\n")
        ends = ".ENDS PORT_{0}\n\n"
        f.write(ends.format(u))
    #END MAIN SUBCIRCUIT - COMMENT IF USING PSPICE
    #f.write(".ENDS\n")
    #Termination
    f.write(".END")
    #Close the file
    f.close()
    
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