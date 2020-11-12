# A simple python script to plot the GW
# signals over time, for a chosen mode

import numpy as np;
import matplotlib.pyplot as plt;

# output data from running merger
M = 1.0
mu = 0.05
v = 0.5
r = 3000
symmetry = 4

# make the plot
fig = plt.figure()

# for v=0.4
data1 = np.loadtxt("Run8/RhoIntegral.dat")
labelstring = "integral Source dV v = 0.8"
timedata = data1[:,0]
Fdata = data1[:,1]
SourceData = Fdata * symmetry
plt.plot(timedata, SourceData, '--', lw = 1.0, label=labelstring)

# for v=0.4
data1 = np.loadtxt("Run7/RhoIntegral.dat")
labelstring = "integral Source dV v = 0.7"
timedata = data1[:,0]
Fdata = data1[:,1]
SourceData = Fdata * symmetry
plt.plot(timedata, SourceData, '--', lw = 1.0, label=labelstring)

# for v=0.4
data1 = np.loadtxt("Run6/RhoIntegral.dat")
labelstring = "integral Source dV v = 0.6"
timedata = data1[:,0]
Fdata = data1[:,1]
SourceData = Fdata * symmetry
plt.plot(timedata, SourceData, '--', lw = 1.0, label=labelstring)

# mass dataset - source
data1 = np.loadtxt("Run5/RhoIntegral.dat")
labelstring = "integral Source dV v = 0.5"
timedata = data1[:,0]
Fdata = data1[:,1]
SourceData = Fdata * symmetry
plt.plot(timedata, SourceData, '--', lw = 1.0, label=labelstring)

# for v=0.4
data1 = np.loadtxt("Run4/RhoIntegral.dat")
labelstring = "integral Source dV v = 0.4"
timedata = data1[:,0]
Fdata = data1[:,1]
SourceData = Fdata * symmetry
plt.plot(timedata, SourceData, '--', lw = 1.0, label=labelstring)

# for v=0.3
data1 = np.loadtxt("Run3/RhoIntegral.dat")
labelstring = "integral Source dV v = 0.3"
timedata = data1[:,0]
Fdata = data1[:,1]
SourceData = Fdata * symmetry
plt.plot(timedata, SourceData, '--', lw = 1.0, label=labelstring)

# make the plot look nice
plt.xlabel("time")
plt.ylabel("F")
#plt.xlim(0, 1000)
#plt.ylim(1e-1, 1e2)
plt.legend(loc=3)

# save as png image
filename = "FvsT" + "_mu" + str(mu) + ".png"
plt.savefig(filename)