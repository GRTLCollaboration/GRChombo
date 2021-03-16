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
# Source
Source = symmetry*data1[:,1]
#xMom
dxMom = symmetry*data1[:,2] - symmetry*data1[0,2]

# flux dataset out
data1 = np.loadtxt("Dina/SurfaceIntegrals.dat")
labelstring = "integral(Flux * dt)"
timedata = data1[:,0]
dt = timedata[1] - timedata[0]
NetMomoFlux = data1[:,3]
NetMomiFlux = data1[:,1]
FMomodt = np.zeros_like(timedata)
FMomidt = np.zeros_like(timedata)
Source_dt = np.zeros_like(timedata)
for i, t in enumerate(timedata) :
    if (i > 0) : #and (i < np.size(timedata) - 1) :
        FMomodt[i] += FMomodt[i-1] + NetMomoFlux[i] * dt
        FMomidt[i] += FMomidt[i-1] + NetMomiFlux[i] * dt
        Source_dt[i] += Source_dt[i-1]+ Source[i] * dt

plt.plot(timedata, dxMom, '-', lw = 1.0, label="d(xMom)")
plt.plot(timedata, Source_dt, '-', lw = 1.0, label="Source dt")
plt.plot(timedata, FMomodt, '-', lw = 1.0, label="Flux o dt")
plt.plot(timedata, FMomidt, '-', lw = 1.0, label="Flux i dt")
plt.plot(timedata, FMomodt - FMomidt - Source_dt, '--', lw = 1.0, label="check (flux - source) = dxMom")


# make the plot look nice
plt.xlabel("time")
plt.ylabel("Change in Cloud xMom")
#plt.xlim(0, 100)
#plt.ylim(-10, 10)
plt.legend(loc=0)
plt.grid()

# save as png image
filename = "MvsT" + "_mu" + str(mu) + ".png"
plt.savefig(filename)
