#!/usr/bin/env python
"""
Created on Mon 2 Dec 2013

Script viser import av funksjoner fra numpy og bruk av noen.

@author Benedicte Emilie Braekken
"""
from numpy import *

print 'e^1 =', exp( 1 )                     # Eksponentialfunksjonen
print 'cos(pi) =', cos( pi )                # Cosinus
print 'sqrt(4) =', sqrt( 4 )                # Kvadratrot
print 'range(5) =', range(5)                # Rekke opp til 4
print 'zeros(5) =', zeros(5)                # Tom array med 5 elementer
print 'linspace(0,5,5) =', linspace(0,5,5)  # Rekke som ikke oeker med 1

"""
bruker @ unix $ python numpy_functions.py
e^1 = 2.71828182846
cos(pi) = -1.0
sqrt(4) = 2.0
range(5) = [0, 1, 2, 3, 4]
zeros(5) = [ 0.  0.  0.  0.  0.]
linspace(0,5,5) = [ 0.    1.25  2.5   3.75  5.  ]
"""
