# A simple python script to plot the GW
# signals over time, for a chosen mode

import numpy as np;
import matplotlib.pyplot as plt;

# output data from running merger
M = 1.0
mu = 0.05
v = 0.0
r = 500
data1 = np.loadtxt("RhoIntegral.dat")
data2 = np.loadtxt("Force_integrals.dat")

# make the plot
fig = plt.figure()

# first dataset
labelstring = "M-M0, v = " + str(v) + " mu = " + str(mu) + " r = " + str(r)
timedata = data1[:,0]
Mdata = data1[:,2] - data1[0,2]
#plt.plot(timedata, Mdata, '-', lw = 1.0, label=labelstring)

# second dataset net
labelstring = "integral net Flux dt$"
timedata = data2[:,0]
Fdata = (data2[:,3] - data2[:,1])
deltaE = np.zeros_like(Fdata)
for i, F in enumerate(Fdata) :
   if i==0 :
      deltaE[i] = 0.0
   else :
      deltaE[i] = deltaE[i-1] + F * (timedata[i] - timedata[i-1])
plt.plot(timedata, deltaE, '-', lw = 1.0, label=labelstring)

# second dataset out
labelstring = "integral out Flux dt$"
timedata = data2[:,0]
Fdata = data2[:,3]
deltaE = np.zeros_like(Fdata)
for i, F in enumerate(Fdata) :
   if i==0 :
      deltaE[i] = 0.0
   else :
      deltaE[i] = deltaE[i-1] + F * (timedata[i] - timedata[i-1])
plt.plot(timedata, deltaE, '-', lw = 1.0, label=labelstring)

# second dataset in
labelstring = "integral in Flux dt$"
timedata = data2[:,0]
Fdata = data2[:,1]
deltaE = np.zeros_like(Fdata)
for i, F in enumerate(Fdata) :
   if i==0 :
      deltaE[i] = 0.0
   else :
      deltaE[i] = deltaE[i-1] + F * (timedata[i] - timedata[i-1])
plt.plot(timedata, deltaE, '-', lw = 1.0, label=labelstring)

# third dataset
labelstring = "integral Source + DeltaM dt$"
timedata = data1[:,0]
Fdata = data1[:,1]
deltaE = np.zeros_like(Fdata)
for i, F in enumerate(Fdata) :
   if i==0 :
      deltaE[i] = 0.0
   else :
      deltaE[i] = deltaE[i-1] + F * (timedata[i] - timedata[i-1])
plt.plot(timedata, deltaE + Mdata, '-', lw = 1.0, label=labelstring)

# make the plot look nice
plt.xlabel("time")
plt.ylabel("xMom")
#plt.xlim(0, 1000)
#plt.ylim(1e-1, 1e2)
plt.legend()

# save as png image
filename = "EvsT" + "_mu" + str(mu) + ".png"
plt.savefig(filename)