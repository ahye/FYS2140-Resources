#!/usr/bin/env python
"""
Created on Mon 2 Dec 2013

Scriptet viser hvordan man kan plotte 3D-figurer med matplotlib.

@author Benedicte Emilie Braekken
"""
from matplotlib.pyplot import *
from numpy import *
from mpl_toolkits.mplot3d import Axes3D

def gauss_2d( x, y ):
    '''
    Todimensjonal gauss-kurve.
    '''
    # Konstant
    A = 1

    X = ( x - x0 )**2 / ( 2. * sigma_x**2 )
    Y = ( y - y0 )**2 / ( 2. * sigma_y**2 )

    return A * exp( - ( X + Y ) )

# Antall punkter hver retning (reelt blir det n^2 for hele rommet)
n = 3e2

# Senteret gauss-kurven ligger paa
x0 = 0
y0 = 0

# Bredden paa gauss-kurven
sigma_x = 1
sigma_y = 1

# Enhetsvektorne i hver retning
x = linspace( x0 - 4. * sigma_x, x0 + 4. * sigma_x, n )
y = linspace( y0 - 4. * sigma_y, y0 + 4. * sigma_y, n )

# Lager de to tabellene som inneholder hvert punkt i rommet
# de to enhetsvektorne over utspenner
X, Y = meshgrid( x, y )

# Lager figur
fig = figure()

# Lager akser
ax = fig.add_subplot( 111, projection='3d' )

# Lager overflateplott
ax.plot_surface( X, Y, gauss_2d( X, Y ) )

# Viser
show()
