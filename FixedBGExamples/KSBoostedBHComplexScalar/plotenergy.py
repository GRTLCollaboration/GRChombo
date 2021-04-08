# A simple python script to plot the GW
# signals over time, for a chosen mode

import numpy as np;
import matplotlib.pyplot as plt;

# output data for setup
M = 1.0
mu = 0.5
r = 300
symmetry = 4
# make the plot
fig = plt.figure()

# volume integral dataset out
data1 = np.loadtxt("Dina/RhoIntegral.dat")
timedata = data1[:,0]
dM = symmetry*data1[:,3] - symmetry*data1[0,3]
#print(dM)

# flux dataset out
data1 = np.loadtxt("Dina/SurfaceIntegrals.dat")
labelstring = "integral(Flux * dt)"
timedata = data1[:,0]
dt = timedata[1] - timedata[0]
NetEoFlux = data1[:,4]
NetEiFlux = data1[:,2]
FEodt = np.zeros_like(timedata)
FEidt = np.zeros_like(timedata)
for i, F in enumerate(timedata) :
    if (i > 0) :
        FEodt[i] += FEodt[i-1] + NetEoFlux[i] * dt
        FEidt[i] += FEidt[i-1] + NetEiFlux[i] * dt

plt.plot(timedata, FEodt, '-', lw = 1.0, label="Edot outer dt")
plt.plot(timedata, FEidt, '-', lw = 1.0, label="Edot inner dt")
plt.plot(timedata, dM, '-', lw = 1.0, label="M-M0")
plt.plot(timedata, FEodt - FEidt, '--', lw = 1.0, label="check M-M0")

# make the plot look nice
plt.xlabel("time")
plt.ylabel("Change in Cloud Mass")
#plt.xlim(0, 100)
#plt.ylim(-10, 10)
plt.legend(loc=0)
plt.grid()

# save as png image
filename = "EvsT" + "_mu" + str(mu) + ".png"
plt.savefig(filename)
