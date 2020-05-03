=== Vector Fitting Algorithm in Python ===

Contributor: Jared Reichle
Requires: Python 3.0
Tested up to: Python 3.7
Created: January 23rd, 2020

=== Description ===

The Vector Fitting Algoirthm in Python is a python script meant
to take in fequency data and return a transfer function. That
function is meant to show the magnitude/phase relationship that
reflects the overall network.

In this program, there are different steps, or functions, on 
which the program runs. The functions include the main function,
vectfit_step function, calculate_residues function, and the
vectfit_auto function.

The program will take an 18th order system with varrying
frequency data points (both complex and real) and will find the
best fit for the data. The deliverble should be near matched with
the given poles-residues and the input transfer function.

=== Instructions ===

Open ex3.py or ex4.py along with Vectorfitting.py. You can run
ex1.py or ex2.py with the given parameters or change the parameters
to a desired system. ex3.py runs a single port system while ex4.py
runs a 4 port system.

=== Parameters ===

s -- Frequency samples: These are meant to be random points of data
that represent the frequency response of the system.
f -- Frequency response f(s): This is meant to be the model system
that our vector fitting will be fit off of.
d -- Offset: This is the starting magnitude value
h -- Slope: This is the slope of the fitted vector and how it fits
with the whole system.


=== Contact ===

If you have any questions or suggestions, you can contact me
below:
Work email: reic2736@vandals.uidaho.edu
Personal email: reichle.jared@gmail.com