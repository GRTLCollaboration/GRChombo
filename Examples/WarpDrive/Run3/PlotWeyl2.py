# A simple python script to plot the GW
# signals over time, for a chosen mode

import numpy as np;
import matplotlib.pyplot as plt;
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import sys

mode = str(sys.argv[1])
rad = int(sys.argv[2])

# coord locations of extraction radii
R0 = 2.0 #200
R1 = 2.0 #300
# Total ADM mass
M = 1.0
# The mode, as text
#mode = "21"
# output data from running merger
data = np.loadtxt("Weyl_integral_" + mode + ".dat")
data2 = np.loadtxt("../Run2/Weyl_integral_" + mode + ".dat")

# make the plot
fig = plt.figure()

# second radius
r1 = 100 #R1 + M*np.log(R1/(2.0*M) - 1.0)
timedata1 = (data[:,0] - r1) / M
fluxdata1 = data[:,rad*2-1]
plt.plot(timedata1, fluxdata1, '-', lw = 0.75, label="HR_r"+str(rad))

# third radius
r1 = 100 #R1 + M*np.log(R1/(2.0*M) - 1.0)
timedata1 = (data2[:,0] - r1) / M
fluxdata1 = data2[:,rad*2-1]
plt.plot(timedata1, fluxdata1, '--', lw = 0.75, label="LR_r"+str(rad))

# make the plot look nice
plt.xlabel("time t / M")
plt.ylabel("Re (Psi4) el, em = " + mode)
plt.xlim(-20, 600)
#plt.ylim(1.85, 1.88)
plt.legend()

# save as png image
filename = "Weyl_" + mode + ".png"
plt.savefig(filename)
