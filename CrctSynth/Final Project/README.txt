=== FINAL PROJECT ===

Contributor: Jared Reichle
Requires: Python 3.0
Tested up to: Python 3.7
Created: May 1st, 2020

=== Description ===

This program utilizes the fast-relaxed vector fitting algorithm and the
passivity enforcement algorithm as outlined in the references laid out in
Final_Project.py. The program will take s-parameter data and bigY data from
ex2_Y.mat file and then create a state-space and pole-residue model that best
fits the data with an N-order model. That N order model is then checked for
passivity violations, perturbated, and passivity is enforced for simulation
stability. The stable, passive system is then written to a spice netlist for
simulation with Advance Deisgn Systems (ADS).

=== Instructions ===

Ensure that the following files are all within your working directory:
    -Final_Project.py
    -VFDriver.py
    -RPDriver.py
    -vectfit3.py
    -intergheig.py
    -fitcalcABCDE.py
    -fitcalcPRE.py
    -rot.py
    -pr2ss.py
    -FRPY.py
    -violextrema.py
    -ex2_Y.mat
    -README.txt
    
Then, open up Final_Project.py and execute program. This should show the
results from the vector fitting and passivity enforcement in real time, and
produce a SPICE netlist within the working directory labeled spice_net_list.sp

=== Parameters ===

bigY - The admittance matrix size Np x Np x D where Np is the number of ports
in the system and D is the number of frequency data points.

s - The row vector of the frequency data points in a Laplace form (s = j(omega))

The state space model:

          N   SERCijm
  Yij(s)=SUM(---------) +SERDij +s*SEREij ("ij" denotes element i,j)
         m=1 (s-SERAm)

=== Returns ===

spice_net_list.sp - The spice netlist that has been passively enforced and fitted to the data

bigYfit - The fitted data to bigY. Accuracy is dependent on number of iterations
and order (N) of model.

rmserr - The root mean square error of the whole system to see how close your
fitted and perturbated model comes to the original data

=== Contact ===

If you have any questions or suggestions, you can contact me
below:
Work email: reic2736@vandals.uidaho.edu
Personal email: reichle.jared@gmail.com