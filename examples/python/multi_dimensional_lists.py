#!/usr/bin/env python
"""
Created on Mon 2 Dec 2013

Kort eksempel som viser flerdimensjonale lister og arrayer.

@author Benedicte Emilie Braekken
"""
from numpy import *

# Ren python liste
min_liste = [ [1, 2, 3], [4, 5, 6], [7, 8, 9] ]
print min_liste

# Indeksering
print 'Matrise som liste:', min_liste[0]
print 'Topp venstre element:', min_liste[0][0]

# Numpy array
min_array = array(min_liste)
print min_array

# Indeksering
print 'Matrise som numpy-array:', min_array[0]
print 'Topp venstre element:', min_array[0,0]

"""
bruker @ unix $ python multi_dimensional_list.py
[[1, 2, 3], [4, 5, 6], [7, 8, 9]]
Matrise som liste: [1, 2, 3]
Topp venstre element: 1
[[1 2 3]
 [4 5 6]
 [7 8 9]]
Matrise som numpy-array: [1 2 3]
Topp venstre element: 1
"""
