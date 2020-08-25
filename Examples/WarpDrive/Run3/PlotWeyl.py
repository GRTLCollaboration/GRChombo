# A simple python script to plot the GW
# signals over time, for a chosen mode

import numpy as np;
import matplotlib.pyplot as plt;
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import sys

mode = str(sys.argv[1])

# coord locations of extraction radii
R0 = 2.0 #200
R1 = 2.0 #300
# Total ADM mass
M = 1.0
# The mode, as text
#mode = "21"
# output data from running merger
data = np.loadtxt("Weyl_integral_" + mode + ".dat")

# make the plot
fig = plt.figure()

# first radius
r0 = 0 #R0 + M*np.log(R0/(2.0*M) - 1.0)
timedata0 = (data[:,0] - r0) / M
fluxdata0 = data[:,1]
plt.plot(timedata0, fluxdata0, ':', lw = 0.75, label="r0")

# second radius
r1 = 100 #R1 + M*np.log(R1/(2.0*M) - 1.0)
timedata1 = (data[:,0] - r1) / M
fluxdata1 = data[:,3]
plt.plot(timedata1, fluxdata1, '--', lw = 0.75, label="r1")

# third radius
r2 = 200 #R1 + M*np.log(R1/(2.0*M) - 1.0)
timedata2 = (data[:,0] - r2) / M
fluxdata2 = data[:,5]
plt.plot(timedata2, fluxdata2, '-', lw = 0.75, label="r2")

# 4th radius
r3 = 300 #R1 + M*np.log(R1/(2.0*M) - 1.0)
timedata3 = (data[:,0] - r3) / M
fluxdata3 = data[:,7]
plt.plot(timedata3, fluxdata3, '-', lw = 0.75, label="r3")

# 5th radius
r4 = 400 #R1 + M*np.log(R1/(2.0*M) - 1.0)
timedata4 = (data[:,0] - r2) / M
fluxdata4 = data[:,9]
plt.plot(timedata4, fluxdata4, '-', lw = 0.75, label="r4")

# make the plot look nice
plt.xlabel("time t / M")
plt.ylabel("Re (Psi4) el, em = " + mode)
plt.xlim(-20, 600)
#plt.ylim(1.85, 1.88)
plt.legend()

# save as png image
filename = "Weyl_" + mode + ".png"
plt.savefig(filename)
