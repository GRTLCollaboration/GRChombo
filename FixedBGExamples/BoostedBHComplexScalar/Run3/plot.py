# A simple python script to plot the GW
# signals over time, for a chosen mode

import numpy as np;
import matplotlib.pyplot as plt;

# output data from running merger
M = 1.0
mu = 0.05
v = 0.10
labelstring = "v=" + str(v) + " mu=" + str(mu)
data1 = np.loadtxt("Force_integrals.dat")

# make the plot
fig = plt.figure()

# first dataset
timedata = data1[:,0] / M
Fdata = data1[:,1]
plt.plot(timedata, Fdata, '-', lw = 1.0, label=labelstring+"r=100")

Fdata = data1[:,3]
plt.plot(timedata, Fdata, '-', lw = 1.0, label=labelstring+"r=200")

Fdata = data1[:,5]
plt.plot(timedata, Fdata, '-', lw = 1.0, label=labelstring+"r=300")

# make the plot look nice
plt.xlabel("time t [M]")
plt.ylabel("F")
#plt.xlim(0, 1000)
#plt.ylim(1e-1, 1e2)
plt.legend(loc=0)

# save as png image
filename = "FvsT_v" + str(v) + "_mu" + str(mu) + ".png"
plt.savefig(filename)