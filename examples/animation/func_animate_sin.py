#!/usr/bin/env python
"""
Created on Mon 2 Dec 2013

Eksempelscript som viser hvordan en sinusboelge kan animeres med
funksjonsanimasjon.

@author Benedicte Emilie Braekken
"""
from numpy import *
from matplotlib.pyplot import *
from matplotlib import animation

def wave( x, t ):
    '''
    Funksjonen beskriver en sinusboelge ved tiden t og punktet x.
    '''
    omega = 1   # Vinkelhastighet
    k = 1       # Boelgetall

    return sin( k * x - omega * t )

T = 10
dt = 0.01
nx = 1e3
nt = int( T / dt ) # Antall tidssteg
t = 0

all_waves = [] # Tom liste for aa ta vare paa boelgetilstandene
x = linspace( -pi, pi, nx )

while t < T:
    # Legger til en ny boelgetilstand for hver kjoering
    all_waves.append( wave( x, t ) )

    t += dt

# Tegner initialtilstanden
fig = figure() # Passer paa aa ta vare paa figuren
line, = plot( x, all_waves[0] )
draw()

# Konstanter til animasjonen
FPS = 60 # Bilder i sekundet
inter = 1. / FPS # Tid mellom hvert bilde

def init():
    '''
    '''
    line.set_data( [], [] )
    return line,

def get_frame( frame ):
    '''
    '''
    line.set_data( x, all_waves[ frame ] )
    return line,

anim = animation.FuncAnimation( fig, get_frame, init_func=init,
                                frames=nt, interval=inter, blit=True )

show()
