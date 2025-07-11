=== Passivity Enforcement in Python ===

Contributor: Jared Reichle
Requires: Python 3.0
Tested up to: Python 3.7
Created: April 12th, 2020

=== Description ===

This program exemplifies the implementation of the passivity
enforcement algorithm of a state-space model for an arbitrary
N-port network. This specific program does the implementation for
an N-port network using the Y parameter (Admittance matrix)

In this example, I include the A, C, D, and E matrices that make
up the vector-fitted (Project 1) state-space model using varying
frequency data. This state space model is in the form of .csv files
within the working directory. This program will then plot the
eigenvalues and print off a SPICE netlist that should run stable
within the time domain simulation.

Testing for the time domain stability and the assurance of passivity
can be done through ADS. This program does not implement the testing.

=== Instructions ===

Ensure that the following files are all within your working directory:
    -SERA.csv
    -SERC.csv
    -SERD.csv
    -SERE.csv
    -sPR.csv
    -Project3.py
    -README.txt
    
Then, open up Project3.py and execute program. This should show the
results from the passivity enforcement algorithm and product a SPICE
netlist titled P3.sp. This can then be imported into ADS and tested.

=== Parameters ===

SERA - A portion of the state space model
SERC - C portion of the state space model
SERD - D portion of the state space model
SERE - E portion of the state space model
sPR - frequency samples matrix
bigY - The large Y-parameter network

The state space model:

          N   SERCijm
  Yij(s)=SUM(---------) +SERDij +s*SEREij ("ij" denotes element i,j)
         m=1 (s-SERAm)

=== Contact ===

If you have any questions or suggestions, you can contact me
below:
Work email: reic2736@vandals.uidaho.edu
Personal email: reichle.jared@gmail.com