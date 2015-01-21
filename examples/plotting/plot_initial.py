#!/usr/bin/env python
"""
Created on Mon 2 Dec 2013

Plotter grunntilstanden for en harmonisk oscillator.

@author Benedicte Emilie Brakken
"""
from numpy import *
from matplotlib.pyplot import *

# Kun for aa teste
omega = 1           # [rad / s]

# Fysiske parametre
hbarc = 0.1973      # [MeV pm]
E0p = 938.27        # [MeV]
c = 3e2             # [pm / as]

# x-verdier
x = linspace( -pi, pi, 1e4 )

def Psi0( x ):
    '''
    Grunntilstanden for en harmonisk oscillator.
    '''
    A = ( E0p * omega / ( pi * hbarc * c ) )**0.25
    B = exp( - E0p * omega / ( 2 * hbarc * c) * x**2 )
    return A * B

# Henter funksjonsverdier og lagrer i arrayen Psi
Psi = Psi0(x)

# Lager et nytt figurvindu
figure()

# Plotter x mot Psi0
plot( x, abs( Psi )**2 )

# Tekst langs x-aksen
xlabel('$x$ [pm]')

# Tekst langs y-aksen
ylabel('$|\Psi_0 (x, 0)|^2$ [1/pm]')

# Tittel paa plottet
title('Grunntilstanden for harmonisk oscillator')

# Viser det vi har plottet
show()
