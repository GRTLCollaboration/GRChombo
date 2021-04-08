# A simple python script to plot the GW
# signals over time, for a chosen mode

import numpy as np;
import matplotlib.pyplot as plt;

# output data for setup
M = 1.0
mu = 0.05
r = 300
symmetry = 4
# make the plot
fig = plt.figure()

# volume integral dataset out
data1 = np.loadtxt("VolumeIntegrals.dat")
timedata = data1[:,0]
dM = symmetry*data1[:,1] - symmetry*data1[0,1]
Source = symmetry*data1[:,2]


#data1 = np.loadtxt("VolumeIntegrals.dat")
#dM2 = symmetry*data1[:,1] - symmetry*data1[0,1]

# flux dataset out
data1 = np.loadtxt("SurfaceIntegrals.dat")
labelstring = "integral(Flux * dt)"
timedata = data1[:,0]
dt = timedata[1] - timedata[0]
NetEoFlux = data1[:,4]
NetEiFlux = data1[:,1]
FEodt = np.zeros_like(timedata)
FEidt = np.zeros_like(timedata)
Source_dt = np.zeros_like(timedata)
for i, F in enumerate(timedata) :
    if (i > 0) :
       FEodt[i] += FEodt[i-1] + NetEoFlux[i] * dt
       FEidt[i] += FEidt[i-1] + NetEiFlux[i] * dt
       Source_dt[i] += Source_dt[i-1]+ Source[i] * dt

plt.plot(timedata, FEodt, '-', lw = 1.0, label="Edot outer dt")
plt.plot(timedata, FEidt, '-', lw = 1.0, label="Edot inner dt")
plt.plot(timedata, Source_dt, '-', lw = 1.0, label="Source")
plt.plot(timedata, dM, '-', lw = 1.0, label="M-M0")
#plt.plot(timedata, dM2, '-', lw = 1.0, label="M-M0 LL")
plt.plot(timedata, FEodt - FEidt + Source_dt, '--', lw = 1.0, label="check M-M0")

# make the plot look nice
plt.xlabel("time")
plt.ylabel("Change in Cloud Mass")
#plt.xlim(0, 100)
#plt.ylim(-10, 10)
plt.legend(loc=0)
plt.grid()

# save as png image
filename = "EvsT.png"
plt.savefig(filename)
