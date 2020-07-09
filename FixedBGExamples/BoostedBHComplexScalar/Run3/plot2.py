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
#N = 36 #width of rolling average


# make the plot
fig = plt.figure()

# first dataset
timedata = data1[:,0] / M
Fdata = data1[:,1]

length = np.size(timedata)
#total = 0.0
#count = 0
#time = np.zeros(length-N)
#Fmean = np.zeros(length-N)
#for i, t in enumerate(timedata) :
#    if (i < N/2) :
#        total += Fdata[i]
#    elif (length - i - 1) < N/2 :
#        break
#    else :
#        total += Fdata[i+N/2] - Fdata[i-N/2]
#        time[count] = timedata[i]
#        Fmean[count] = total/N
#        count += 1
#        if (Fdata[i] * Fdata[i-1] < 0) :
#            print ("Switch at datapt ", i, " time ", timedata[i])


integralF = np.zeros_like(timedata)
for i, t in enumerate(timedata) :
    if (i ==0) :
        integralF[i] = Fdata[i]
    else :
        integralF[i] = integralF[i-1] + Fdata[i]

#plt.plot(timedata, integralF, '-', lw = 1.0, label=labelstring+"r=200")

Fdata = data1[:,3]
integralF = np.zeros_like(timedata)
for i, t in enumerate(timedata) :
    if (i ==0) :
        integralF[i] = Fdata[i]
    else :
        integralF[i] = integralF[i-1] + Fdata[i]
plt.plot(timedata, integralF, '-', lw = 1.0, label=labelstring+"r=300")

Fdata = data1[:,5]
integralF = np.zeros_like(timedata)
for i, t in enumerate(timedata) :
    if (i ==0) :
        integralF[i] = Fdata[i]
    else :
        integralF[i] = integralF[i-1] + Fdata[i]
#plt.plot(timedata, integralF, '-', lw = 1.0, label=labelstring+"r=400")

# make the plot look nice
plt.xlabel("time t [M]")
plt.ylabel("Fmean")
#plt.xlim(0, 1000)
#plt.ylim(-1, 1)
plt.legend(loc=0)

# save as png image
filename = "FmeanvsT_v" + str(v) + "_mu" + str(mu) + ".png"
plt.savefig(filename)