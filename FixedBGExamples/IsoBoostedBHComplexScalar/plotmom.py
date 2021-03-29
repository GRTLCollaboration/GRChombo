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
dM = symmetry*data1[:,3] - symmetry*data1[0,3]
Source = symmetry*data1[:,4]

# flux dataset out
data1 = np.loadtxt("SurfaceIntegrals.dat")
labelstring = "integral(Flux * dt)"
timedata = data1[:,0]
dt = timedata[1] - timedata[0]
NetEiFlux = data1[:,3]
NetEoFlux = data1[:,6]
FEodt = np.zeros_like(timedata)
FEidt = np.zeros_like(timedata)
Source_dt = np.zeros_like(timedata)
for i, F in enumerate(timedata) :
    if (i > 0) :
       FEodt[i] += FEodt[i-1] + NetEoFlux[i] * dt
       FEidt[i] += FEidt[i-1] + NetEiFlux[i] * dt
       Source_dt[i] += Source_dt[i-1]+ Source[i] * dt

plt.plot(timedata, FEodt, '-', lw = 1.0, label="Mdot outer dt")
plt.plot(timedata, FEidt, '-', lw = 1.0, label="Mdot inner dt")
plt.plot(timedata, Source_dt, '-', lw = 1.0, label="Source dt")
plt.plot(timedata, dM, '-', lw = 1.0, label="M-M0")
plt.plot(timedata, FEidt - FEodt + Source_dt, '--', lw = 1.0, label="check M-M0")

# make the plot look nice
plt.xlabel("time")
plt.ylabel("Change in Cloud Mom")
#plt.xlim(0, 100)
#plt.ylim(-10, 10)
plt.legend(loc=0)
plt.grid()

# save as png image
filename = "MvsT.png"
plt.savefig(filename)
