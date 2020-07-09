# A simple python script to plot the GW
# signals over time, for a chosen mode

import numpy as np;
import matplotlib.pyplot as plt;

# output data from running merger
M = 1.0
mu = 0.05
data1 = np.loadtxt("Run1/Force_integrals.dat")
data2 = np.loadtxt("Run2/Force_integrals.dat")

# make the plot
fig = plt.figure()
irad = 1
r = 100
Area = 4 * np.pi * r * r

# first dataset
v=0.1
labelstring = "v = " + str(v) + " mu = " + str(mu) + " r = " + str(r)
timedata = data1[:,0] / M
Fdata = data1[:,2*irad - 1]
rhodata = data1[:,2*irad]/Area
plt.plot(timedata, Fdata/rhodata, '-', lw = 1.0, label=labelstring)

# second dataset
v=0.15
labelstring = "v = " + str(v) + " mu = " + str(mu) + " r = " + str(r)
timedata = data2[:,0] / M
Fdata = data2[:,2*irad - 1]
rhodata = data2[:,2*irad]/Area
plt.plot(timedata, Fdata/rhodata, '-', lw = 1.0, label=labelstring)

# make the plot look nice
plt.xlabel("time t [M]")
plt.ylabel("F")
#plt.xlim(0, 1000)
#plt.ylim(1e-1, 1e2)
plt.legend(loc=0)

# save as png image
filename = "FvsT_v" + "_mu" + str(mu) + ".png"
plt.savefig(filename)