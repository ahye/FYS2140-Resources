#!/usr/bin/env python
"""
Created on Sun 27 Oct 2013.

Uses the Crank Nicolson scheme to solve the time dependent Schrodinger equation
for a potential well.

Animation is done using the matplotlib.pyplot library.

Usage:
    python CrankNicolsonPotentialWell.py
or equivalent:
    ./CrankNicolsonPotentialWell.py
No commandline arguments are needed.

@author Benedicte Emilie Braekken
"""
# Tools for sparse matrices
import scipy.sparse as sparse
import scipy.sparse.linalg

# Numerical tools
from numpy import *

# Plotting library
from matplotlib.pyplot import *

"""Physical constants"""
_E0p = 938.27        # Rest energy for a proton [MeV]
_hbarc = 0.1973      # [MeV pm]
_c = 3.0e2           # Spees of light [pm / as]

def Psi0( x ):
    '''
    Initial state for a travelling gaussian wave packet.
    '''
    x0 = 0 # [pm]
    a = 0.0050 # [pm]

    A = ( 1. / ( 2 * pi * a**2 ) )**0.25
    K1 = exp( - ( x - x0 )**2 / ( 4. * a**2 ) )

    return A * K1

def potentialWell( x, depth=35, width=0.15 ):
    """
    Gives the potential for a potential well of depth depth and width width.

    @param depth Gives the depth of the potential well. Given as the magnitude
    (positive integer / double / float).
    @param width Gives the width of the potential well. Given as positive
    integer definig the fraction of the x spectrum to contain the well. For
    example, 1 will mean that the well covers the whole spectrum and 0.5 that
    it covers half of it.
    """
    # Declare new empty array with same length as x
    potential = zeros( len( x ) )

    potential[ (0.5-0.5*width)*len(potential) : (0.5+0.5*width)*len(potential) ] = -depth

    return potential

if __name__ == '__main__':
    nx = 1001 # Number of points in x direction
    dx = 0.001 # Distance between x points [pm]

    # Use zero as center, same amount of points each side
    a = - 0.5 * nx * dx
    b = 0.5 * nx * dx
    x = linspace( a, b, nx )

    # Time parameters
    T = 0.2 # How long to run simulation [as]
    dt = 1e-4 # The time step [as]
    t = 0
    time_steps = int( T / dt ) # Number of time steps

    # Constants - save time by calculating outside of loop
    k1 = - ( 1j * _hbarc * _c) / (2. * _E0p )
    k2 = ( 1j * _c ) / _hbarc

    # Create the initial state Psi
    Psi = Psi0(x)

    # Create the matrix containing central differences. It it used to
    # approximate the second derivative.
    data = ones((3, nx))
    data[1] = -2*data[1]
    diags = [-1,0,1]
    D2 = k1 / dx**2 * sparse.spdiags(data,diags,nx,nx)

    # Identity Matrix
    I = sparse.identity(nx)

    # Create the diagonal matrix containing the potential.
    V_data = potentialWell(x)
    V_diags = [0]
    V = k2 * sparse.spdiags(V_data, V_diags, nx, nx)

    # Put mmatplotlib in interactive mode for animation
    ion()

    # Setup the figure before starting animation
    fig = figure() # Create window
    ax = fig.add_subplot(111) # Add axes
    line, = ax.plot( x, abs(Psi)**2, label='$|\Psi(x,t)|^2$' ) # Fetch the line object

    # Also draw a green line illustrating the potential
    ax.plot( x, V_data, label='$V(x)$' )

    # Add other properties to the plot to make it elegant
    fig.suptitle("Solution of Schrodinger's equation with potential well") # Title of plot
    ax.grid('on') # Square grid lines in plot
    ax.set_xlabel('$x$ [pm]') # X label of axes
    ax.set_ylabel('$|\Psi(x, t)|^2$ [1/pm] and $V(x)$ [MeV]') # Y label of axes
    ax.legend(loc='best') # Adds labels of the lines to the window
    draw() # Draws first window

    # Time loop
    while t < T:
        """
        For each iteration: Solve the system of linear equations:
        (I - k/2*D2) u_new = (I + k/2*D2)*u_old
        """
        # Set the elements of the equation
    	A = (I - dt/2*(D2 + V))
    	b = (I + dt/2. * (D2 + V)) * Psi

        # Calculate the new Psi
    	Psi = sparse.linalg.spsolve(A,b)

        # Update time
    	t += dt

    	# Plot this new state
    	line.set_ydata( abs(Psi)**2 ) # Update the y values of the Psi line
        draw() # Update the plot

    # Turn off interactive mode
    ioff()

    # Add show so that windows do not automatically close
    show()
