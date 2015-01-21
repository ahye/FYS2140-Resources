#!/usr/bin/env python
"""
Created on Tir 29 Oct 2013.

Script uses the forward euler method to solve the Schrodinger equation for a
potential well.

Animation is done using the matplotlib.pyplot library.

Usage:
    python ForwardEulerPotentialWell.py
or equivalent:
    ./ForwardEulerPotentialWell.py
No commandline arguments are needed.

Notice how low the time step needs to be for the method to be stable. Try
raising the timestep to 1e6 and see what happens. Also, compare with the needed
time step for the Crank Nicolson method.

Another thing to note is that when using this method for integration it is
necessary to normalize the wave function for each time step. This is because
the forward euler scheme does not contain the normalization.

@author Benedicte Emilie Braekken
"""
# Tools for sparse matrices
import scipy.sparse as sparse
import scipy.sparse.linalg

# Integration tool needed for normalization
import scipy.integrate

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
    dt = 1e-7 # The time step [as]
    t = 0
    time_steps = int( T / dt ) # Number of time steps

    # Fetch potential
    V = potentialWell(x)

    # Constants - save time by calculating outside of loop
    k1 = - ( 1j * _hbarc * _c) / (2. * _E0p )
    k2 = ( 1j * _c * V ) / _hbarc

    # Create initial state Psi
    Psi = Psi0(x)

    # Initialize empty array for storing second derivative. `complex128`
    # argument necessary for calculations with complex numbers.
    D2Psi = zeros(nx).astype(complex128)

    # Precalculate indexes for finding second order derivative
    IND = arange(1,nx-1)
    INDP = arange(2,nx)
    INDM = arange(0,nx-2)

    # Put mmatplotlib in interactive mode for animation
    ion()

    # Setup the figure before starting animation
    fig = figure() # Create window
    ax = fig.add_subplot(111) # Add axes
    line, = ax.plot( x, abs(Psi)**2, label='$|\Psi(x,t)|^2$' ) # Fetch the line object

    # Also draw a green line illustrating the potential
    ax.plot( x, V, label='$V(x)$' )

    # Add other properties to the plot to make it elegant
    fig.suptitle("Solution of Schrodinger's equation with potential well") # Title of plot
    ax.grid('on') # Square grid lines in plot
    ax.set_xlabel('$x$ [pm]') # X label of axes
    ax.set_ylabel('$|\Psi(x, t)|^2$ [1/pm] and $V(x)$ [MeV]') # Y label of axes
    ax.legend(loc='best') # Adds labels of the lines to the window
    draw() # Draws first window

    # Dont plot all frames. Creates better "framerate" when dealing with
    # methods requiring low time step.
    PLOT_EVERY = 100
    counter = 0 # Counter for drawing frames

    # Time loop
    while t < T:
        # Second order derivative
        D2Psi[IND] = (Psi[INDP]-2.*Psi[IND]+Psi[INDM])/dx**2

        # Wavefunction after timestep
        Psi = Psi + dt*(k1*D2Psi+k2*Psi)

        # Normalize the new wave since this method doesnt keep normalization
        Psi /= scipy.integrate.simps(abs(Psi)**2, dx=dx)

        # Update time
    	t += dt

        counter += 1

        if counter == PLOT_EVERY:
            # Plot this new state
            line.set_ydata( abs(Psi)**2 ) # Update the y values of the Psi line
            draw() # Update the plot
            counter = 0

    # Turn off interactive mode
    ioff()

    # Add show so that windows do not automatically close
    show()
