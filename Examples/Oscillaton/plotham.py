# A simple python script to plot the GW
# signals over time, for a chosen mode

import numpy as np;
import matplotlib.pyplot as plt;

# output data for setup
M = 1.0
mu = 1.0
r = 60
symmetry = 1
# make the plot
fig = plt.figure()

# volume integral dataset out
data1 = np.loadtxt("VolumeIntegrals.dat")
timedata = data1[:,0]
Hamiltonian = symmetry*data1[:,5]
#print(dM)

plt.plot(timedata, Hamiltonian, '-', lw = 1.0, label="Hamiltonian Constraint")

# make the plot look nice
plt.xlabel("time")
plt.ylabel("Change in Hamiltonian Constraint")
#plt.xlim(0, 100)
#plt.ylim(-10, 10)
plt.legend(loc=0)
plt.grid()

# save as png image
filename = "Hamiltonian" + "_mu" + str(mu) + ".png"
plt.tight_layout()
plt.savefig(filename)
