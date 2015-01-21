#!/usr/bin/env python
"""
Created on Mon 2 Dec 2013

Viser hvordan man kan plotte flere ting i forskjellige vinduer med matplotlib.

@author Benedicte Emilie Braekken
"""
from numpy import *
from matplotlib.pyplot import *

# Definerer to tilfeldig valgte funksjoner
def func1( x ):
    return 4*x**2

def func2( x ):
    return exp( -x )

# Lager x-verdier
x = linspace( 0, 10, 1e3 )

# Lager foerst en figur og plotter i den
figure()
plot( x, func1( x ))
title('$4x^2$')     # Setter tittel

# Lager en figur til og plotter i den
figure()
plot( x, func2( x ))
title('$e^{-x}$')   # Setter tittel

# Viser alle plottende
show()
