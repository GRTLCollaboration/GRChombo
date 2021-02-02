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
data1 = np.loadtxt("RhoIntegral.dat")

# make the plot
fig = plt.figure()

# flux dataset out
labelstring = "Source"
timedata = data1[:,0]
Fdata = symmetry*data1[:,1]
plt.plot(timedata, Fdata, '--', lw = 1.0, label=labelstring)

# make the plot look nice
plt.xlabel("time")
plt.ylabel("Force")
#plt.xlim(0, 1000)
#plt.ylim(1e-1, 1e2)
#plt.legend(loc=2)

# save as png image
filename = "FvsT" + "_mu" + str(mu) + ".png"
plt.savefig(filename)