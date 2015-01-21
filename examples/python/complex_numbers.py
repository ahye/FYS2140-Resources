#!/usr/bin/env python
"""
Created on Mon 2 Dec 2013

Kort script som viser hvordan man kan representere komplekse tall.

@author Benedicte Emilie Braekken
"""
from numpy import *

a = 2 + 1j
b = 2j

print '2 + i + 2i =', a + b
print '2 + i - 2i =', a - b
print '2 + i * 2i =', a * b

# Absoluttverdien
print 'abs(2 + i) =', abs(a)
print 'abs(2 + i)^2 =', abs(a)**2

# Kompleks del
print 'imag(2 + i) =', a.imag

# Real del
print 'real(2 + i) =', a.real

# Kompleks numpy array
complex_array = zeros(5, dtype=complex64)
complex_array[0] = a
complex_array[1] = b
print complex_array

"""
bruker @ unix $ python complex_numbers.py
2 + i + 2i = (2+3j)
2 + i - 2i = (2-1j)
2 + i * 2i = (-2+4j)
abs(2 + i) = 2.2360679775
abs(2 + i)^2 = 5.0
imag(2 + i) = 1.0
real(2 + i) = 2.0
[ 2.+1.j  0.+2.j  0.+0.j  0.+0.j  0.+0.j] 
"""
