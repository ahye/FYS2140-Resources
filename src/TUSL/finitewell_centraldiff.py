#####################################################################################
### Program to find eigenenergies of a finite well (dimensionless variables)
### 
### Using the naive approach with the central difference
### to find the next value of the wave function when the potential is defined
### as a finite square well.
###
######################################################################################


# Importing useful stuff
from numpy import *
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
import scipy.integrate
import numpy as np

# Defining finite well potential
def finite_well(W0,z):
	# Well centered around origin of width 2 and depth W0
	zlen = len(z)
	potential = np.zeros(zlen)
	for i in range(0,len(z)):
		if abs(z[i]) < 1:
			# if we are on the domain of the well; make the value -W0
			potential[i] = -W0
	return potential


N = int(1e2)						# number of points
zstart = -1.5
zend = 1.5
z = np.linspace(zstart,zend,N)		# domain of wave function and potential
dz = z[1]-z[0]						# step length
tol = 2.0							# tolerance level 0.01 excludes ground state

Psi = np.zeros(N)					# creating wave function array
Psi[0] = 0							# first point
Psi[1] = 0.1 						# second point (small increase)

hbarc = 197.3 # ev nm
mc2 = 0.511*10**6 # ev

V0 = 77.13
a = 0.2							    # half width of well [nm]
W0 = 2*mc2*a**2*V0/hbarc**2
W = finite_well(W0,z) 		        # getting potential
#V0 = W0*hbarc**2/(2*mc2*a**2)

epsilon_trial = -W0 			# trial eiigenvalue (bound states, need negative energies)
epsilon = []					# list to be filled with epsilon
E = []							# numerical energy list
lastpsi = []					# store last value of Psi
Psi_list = []					# empty list, to be filled with the best Psi

colors = "grcmykw"
color_index = 0
number = 0

# Searching for correct eigenvalue
while epsilon_trial < 0:

	# Calculating all elements in wave function
	for j in range(1,N-1):
		Psi[j+1] = (2 - dz**2*(epsilon_trial - W[j+1]))*Psi[j] - Psi[j-1]
	
	# Normalizing
	Psi /= sqrt(scipy.integrate.simps(abs(Psi)**2,dx=1e-3))

	# Store value of last element in Psi
	Psi_end = abs(Psi[-1])

	# Check if we are within tolerance
	if Psi_end < tol:
		lastpsi.append(Psi_end)
		epsilon.append(epsilon_trial)
		Psi_list.append(list(Psi)) # Add as list to make things behave well

		if len(lastpsi) > 1 and (epsilon[-1] - epsilon[-2]) < 1:
			if lastpsi[-1] < lastpsi[-2]:
				lastpsi.remove(lastpsi[-2])
				epsilon.remove(epsilon[-2])
				Psi_list.remove(Psi_list[-2])

			if lastpsi[-1] > lastpsi[-2]:
				lastpsi.remove(lastpsi[-1])
				epsilon.remove(epsilon[-1])
				Psi_list.remove(Psi_list[-1])

	# Update trial eigenvalue
	epsilon_trial += 0.05


# Physical energies
for i in range(len(epsilon)):
	Energy = epsilon[i]*hbarc**2/(2*mc2*a**2)
	E.append(Energy)


# Print lists of correct eigenvalue and energy
print '-------------------------------------------------------------------------------------------------'
print 'Energy levels of finite potential well of width %.2f nm and depth V0 = -%.2f eV' %(2*a,V0)
print '-------------------------------------------------------------------------------------------------'
print 'Epsilon:' ,epsilon
print 'Numerical energies E [eV]:', E
print '-------------------------------------------------------------------------------------------------'

# Plotting
for i in range(len(Psi_list)):
	Legend = '$\psi_%d$' % (number)
	plot(z,Psi_list[i],color=colors[color_index],label=Legend)
	number += 1
	color_index += 1
plt.title('$Finite \ square \ well$',size=20)
plt.xlabel('$z$',size=18)
plt.ylabel('$\psi(z)$',size=18)
plt.legend(loc='best')
show()