###############################################################################################
## Program finds eigenvalues of the 
## dimensionless harmonic oscillator potential
##
## 				W = z^2
##
## Numerical solutions are plotted with some of the
## known analytic solutions using Hermite polynomials.
## The minus sign in the odd polynomials has no physical meaning.
##
################################################################################################

# Importing useful stuff
from numpy import *
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
import scipy.integrate
import numpy as np

# Defining harmonic oscillator potential
def harmonic_oscillator(z):
	potential = z**2
	return potential

N = int(500)						# number of points
zstart = -5
zend = 5
z = np.linspace(zstart,zend,N)		# domain of wave function and potential
dz = z[1]-z[0]						# step length

Psi = np.zeros(N)
Psi[0] = 0							# first point
Psi[1] = 0.0001 		 			# second point
psi_max = []

epsilon_analytic = [0.5, 1.5, 2.5, 3.5]
epsilon = []						# numerical energy list to be filled with eigenvalues
epsilon_trial = 0.5 			    # trial energy
W = harmonic_oscillator(z) 			# get potential
tol = 5 							# tolerance level (0.01 excludes ground state)

# For plotting numerical solutions with index
colors = 'rcmykwbgrcmykwbgrcmykwbg'
color_index = 0
number = 0

# Searching for correct eigenvalues
while epsilon_trial < 4:
	
	# Calculating next element in wave function
	for j in range(1,N-1):
		Psi[j+1] = (2 - dz**2*(2*epsilon_trial - W[j+1]))*Psi[j] - Psi[j-1]
	
	# Normalizing
	Psi /= sqrt(scipy.integrate.simps(abs(Psi)**2,dx=1e-3))

	if abs(Psi[-1]) < tol:
	 	epsilon.append(epsilon_trial)
		
		# Plot wave function
		figure(1)
		hold('on')
		Label = '$\psi_%d$' % (number)
		plot(z,abs(Psi**2),color=colors[color_index],label=Label)
		number += 1
		color_index += 1

	# Update trial eigenvalue
	epsilon_trial += 0.25


# Print lists of energies
print '-------------------------------------------------------------------------------------------------'
print 'Energy levels of harmonic oscillator potential'
print '-------------------------------------------------------------------------------------------------'
print 'Epsilon numerical: ', epsilon
print 'Epsilon analytical: ', epsilon_analytic
print '-------------------------------------------------------------------------------------------------'


# ANALYTIC SOLUTIONS
# Hermite polynomials
H0 = 1
H1 = 2*z
H2 = 4*z**2 - 2
H3 = 8*z**3 - 12*z

# Symmetric
Psi0 = H0*exp(-z**2/2)
Psi0 /= sqrt(scipy.integrate.simps(abs(Psi0)**2,dx=1e-3))

# Antisymmetric
Psi1 = H1*exp(-z**2/2)
Psi1 /= -sqrt(scipy.integrate.simps(abs(Psi1)**2,dx=1e-3))

# Symmetric
Psi2 = H2*exp(-z**2/2)
Psi2 /= sqrt(scipy.integrate.simps(abs(Psi2)**2,dx=1e-3))

# Antisymmetric
Psi3 = H3*exp(-z**2/2)
Psi3 /= -sqrt(scipy.integrate.simps(abs(Psi3)**2,dx=1e-3))

# PLOTTING
plot(z,abs(Psi0**2),'black',linestyle='--')
plot(z,abs(Psi1**2),'black',linestyle='--')
plot(z,abs(Psi2**2),'black',linestyle='--')
plot(z,abs(Psi3**2),'black',linestyle='--',label='$Analytical$')
plt.ylim(0,12)
plt.xlim(-4,4)
plt.title('$Harmonic \ oscillator$', size=20)
plt.xlabel('$z$', size=18)
plt.ylabel('$|\Psi(z)|^2$', size=18)
legend()
show()