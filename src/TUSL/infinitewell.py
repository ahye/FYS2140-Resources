####################################################################################
###
### Program to find eigenenergies of the infinite square well.
### 
####################################################################################


# Importing useful stuff
from numpy import *
from matplotlib.pyplot import *
import scipy.integrate
import numpy as np
import matplotlib.pyplot as plt

# Defining potential
def infinite_well(z):
	W = zeros(len(z))
	return W

# Constants and parameters
N = 500						# number of points
z = np.linspace(0,1,N)		# position array
dz = z[1]-z[0]				# step length
tol = 0.1					# tolerance level
W = infinite_well(z)		# getting potential
a = 0.4 					# width of well [nm]

hbarc = 197.3 				# eV nm
mc2 = 0.511*10**6 			# eV

Psi = np.zeros(N)			# wave function
Psi[0] = 0 					# initial condition (function must die in endpoints)
Psi[1] = 0.1 				# initial condition
epsilon = [] 				# list to be filled with epsilon
epsilon_anal = []			# analtyic energy list to be filled
E_n = []					# analytical energies
E = []						# numerical energies
lastpsi = []				# value of last psi
Psi_list = []				# list to store the best Psi
epsilon_trial = 9			# trial eigenvalue


# For plotting numerical solutions with index
number = 0 					# in use when labelling wavefunctions in plot
colors = 'cmygbcmygb'	# for different colors in plot
color_index = 0

# Search for correct eigenvalue
while epsilon_trial < 160:
	
	# Calculating wave function
	for j in range(1,N-1):
		Psi[j+1] = (2 - dz**2*(epsilon_trial-W[j+1]))*Psi[j] - Psi[j-1]

	# Normalizing
	Psi /= sqrt(scipy.integrate.simps(abs(Psi)**2,dx=1e-3))

	# Store value of last element in Psi
	Psi_end = abs(Psi[-1])

	# Check if last element is within tolerance
	if Psi_end < tol:
		epsilon.append(epsilon_trial)
		lastpsi.append(Psi_end)
		Psi_list.append(list(Psi)) # add as list to make it behave well

		# Only keep those epsilon and Psi giving minimal value of Psi[-1]
		if len(lastpsi) > 1 and (epsilon[-1] - epsilon[-2]) < 2:
			if lastpsi[-1] < lastpsi[-2]:
				lastpsi.remove(lastpsi[-2])
				epsilon.remove(epsilon[-2])
				Psi_list.remove(Psi_list[-2])
			if lastpsi[-1] > lastpsi[-2]:
				lastpsi.remove(lastpsi[-1])
				epsilon.remove(epsilon[-1])
				Psi_list.remove(Psi_list[-1])

	# Update trial eigenvalue
	epsilon_trial += 0.4

# Physical energies
for i in range(0,len(epsilon)):
	eps = epsilon[i]
	E_phys = eps*hbarc**2/(2*mc2*a**2)
	E.append(E_phys)


# ANALYTIC SOLUTIONS
num = [1,2,3,4]

# Determining energy and wavefunction:
for n in num:
	E_physical = n**2*hbarc**2*pi**2/(2*mc2*a**2)
	E_n.append(E_physical)
	Psi_anal = sin(pi*z*n)

	# Normalizing:
	Psi_anal /= sqrt(scipy.integrate.simps(abs(Psi_anal)**2,dx=1e-3))
	plot(z,Psi_anal,'k--')

# Print lists of energies
print '-------------------------------------------------------------------------------------------------'
print 'Energy levels of infinite potential well of width %.2f nm:' %a
print '-------------------------------------------------------------------------------------------------'
print 'Epsilon: ',epsilon
print 'Numerical energies E [eV]: ', E
print 'Analytical energies En [eV]: ', E_n
print '-------------------------------------------------------------------------------------------------'


# Plotting
for i in range(len(Psi_list)):
	Legend = '$\psi_%d$' % (number)
	plot(z,Psi_list[i],color=colors[color_index],label=Legend)
	number += 1
	color_index += 1

# Axes and title
plt.title('$Infinite \ well$',size=20)
plt.xlabel('$z = x/a$',size=18)
plt.ylabel('$\psi(z)$',size=18)
plot([0,0],[0,0],'k--',label='$Analytical$')
plt.legend(loc='best')
show()