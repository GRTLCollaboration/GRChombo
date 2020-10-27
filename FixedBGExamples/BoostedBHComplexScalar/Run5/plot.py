# A simple python script to plot the GW
# signals over time, for a chosen mode

import numpy as np;
import matplotlib.pyplot as plt;

# output data from running merger
M = 1.0
mu = 0.05
v = 0.1
r = 500
symmetry = 4
data1 = np.loadtxt("RhoIntegral.dat")
data2 = np.loadtxt("Force_integrals.dat")

# make the plot
fig = plt.figure()

# flux dataset out
labelstring = "integral outer Flux dt"
timedata = data2[:,0]
Fdata = data2[:,3]
deltaE = np.zeros_like(Fdata)
for i, F in enumerate(Fdata) :
   if i==0 :
      deltaE[i] = 0.0
   else :
      deltaE[i] = deltaE[i-1] + F * (timedata[i] - timedata[i-1])
plt.plot(timedata, deltaE, '--', lw = 1.0, label=labelstring)
OuterFluxData = deltaE

# mass dataset
labelstring = "M-M0, v = " + str(v) + " mu = " + str(mu) + " r = " + str(r)
timedata = data1[:,0]
Mdata = data1[:,2] - data1[0,2]
Mdata = symmetry * Mdata
plt.plot(timedata, Mdata, '--', lw = 1.0, label=labelstring)

# flux dataset in
labelstring = "integral inner Flux dt"
timedata = data2[:,0]
Fdata = data2[:,1]
deltaE = np.zeros_like(Fdata)
for i, F in enumerate(Fdata) :
   if i==0 :
      deltaE[i] = 0.0
   else :
      deltaE[i] = deltaE[i-1] + F * (timedata[i] - timedata[i-1])
plt.plot(timedata, deltaE, '--', lw = 1.0, label=labelstring)
InnerFluxData = deltaE

# mass dataset - source
labelstring = "integral Source dt"
timedata = data1[:,0]
Fdata = data1[:,1]
deltaE = np.zeros_like(Fdata)
for i, F in enumerate(Fdata) :
   if i==0 :
      deltaE[i] = 0.0
   else :
      deltaE[i] = deltaE[i-1] + F * (timedata[i] - timedata[i-1])
SourceData = deltaE * symmetry
plt.plot(timedata, SourceData, '--', lw = 1.0, label=labelstring)

# combine check
labelstring = "Outerflux = DeltaM (+?) Source - InnerFlux"
Fdata = 1.0 * Mdata + 1.0 * SourceData - 1.0 * InnerFluxData
plt.plot(timedata, Fdata, '-', lw = 1.0, label=labelstring)

# make the plot look nice
plt.xlabel("time")
plt.ylabel("xMom")
#plt.xlim(0, 1000)
#plt.ylim(1e-1, 1e2)
plt.legend(loc=3)

# save as png image
filename = "EvsT" + "_mu" + str(mu) + ".png"
plt.savefig(filename)