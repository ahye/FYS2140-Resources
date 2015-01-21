#!/usr/bin/env python
"""
Created on Mon 2 Dec 2013

Kort eksempelscript som viser hvordan en sinusboelge kan bli animert med
matplotlib.

@author Benedicte Emilie Braekken
"""
from numpy import *
from matplotlib.pyplot import *

def wave( x, t ):
    '''
    Funksjonen beskriver en sinusboelge ved tiden t og punktet x.
    '''
    omega = 1   # Vinkelhastighet
    k = 1       # Boelgetall

    return sin( k * x - omega * t )

n = 1e3
t = 0
T = 10
dt = 0.01
x = linspace( -pi, pi, n )

ion()

figure()
line, = plot( x, wave(x, t) ) # Plotter initialtilstand naar t = 0
xlim([-pi, pi]) # Setter x-aksen til aa gaa fra -pi til pi.
draw()

while t < T:
    line.set_ydata( wave( x, t ) )
    draw()

    t += dt

ioff()
show()
