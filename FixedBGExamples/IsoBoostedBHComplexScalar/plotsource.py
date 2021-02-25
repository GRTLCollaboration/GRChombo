# A simple python script to plot the GW
# signals over time, for a chosen mode

import numpy as np;
import matplotlib.pyplot as plt;

# output data from running merger
M = 1.0
mu = 0.5
v = 0.1
r = 3000
symmetry = 4
# make the plot
fig = plt.figure()

# flux dataset out
data1 = np.loadtxt("Run5/RhoIntegral.dat")
labelstring = "Source 5"
timedata = data1[:,0]
Fdata = symmetry*data1[:,1]
plt.plot(timedata, Fdata, '--', lw = 1.0, label=labelstring)

# flux dataset out
data1 = np.loadtxt("Run4/RhoIntegral.dat")
labelstring = "Source 4"
timedata = data1[:,0]
Fdata = symmetry*data1[:,1]
plt.plot(timedata, Fdata, '--', lw = 1.0, label=labelstring)

# flux dataset out
data1 = np.loadtxt("Run3/RhoIntegral.dat")
labelstring = "Source 3"
timedata = data1[:,0]
Fdata = symmetry*data1[:,1]
plt.plot(timedata, Fdata, '--', lw = 1.0, label=labelstring)

# flux dataset out
data1 = np.loadtxt("Run6/RhoIntegral.dat")
labelstring = "Source 6"
timedata = data1[:,0]
Fdata = symmetry*data1[:,1]
plt.plot(timedata, Fdata, '--', lw = 1.0, label=labelstring)

# make the plot look nice
plt.xlabel("time")
plt.ylabel("Force")
#plt.xlim(0, 1000)
#plt.ylim(1e-1, 1e2)
plt.legend(loc=2)
plt.grid()

# save as png image
filename = "FvsT" + "_mu" + str(mu) + ".png"
plt.savefig(filename)