# A simple python script to plot the GW
# signals over time, for a chosen mode

import numpy as np;
import matplotlib.pyplot as plt;

# output data for setup
M = 1.0
mu = 0.5
r = 400
a = 0.99
symmetry = 2
N = 1 # when to normalise to
alpha = M * mu
r_plus = M + np.sqrt(M*M - a*a)
omega_Re = mu * (1.0 - 0.5 * alpha**2 - 0.33 * alpha **3.0)
omega_Im = 0.99 * alpha**7.0 * (a - 2 * r_plus * omega_Re)

# make the plot
fig = plt.figure()

# volume integral dataset out
data1 = np.loadtxt("Run0/ProcaDensities.dat")
labelstring = "M/M0"
timedata = data1[:,0]
dM = symmetry*data1[:,1]/symmetry*data1[N,1]
plt.semilogy(timedata, dM, '-', lw = 1.0, label=labelstring)

# volume integral dataset out
data1 = np.loadtxt("Run0/ProcaDensities.dat")
labelstring = "J/J0"
timedata = data1[:,0]
dJ = symmetry*data1[:,2]/symmetry*data1[N,2]
plt.semilogy(timedata, dJ, '-', lw = 1.0, label=labelstring)

# analytic
plt.semilogy(timedata, dM[N]*np.exp(2*omega_Im * timedata), '--', lw = 1.0, label="analytic")

# make the plot look nice
plt.xlabel("time")
plt.ylabel("Cloud E/J")
#plt.xlim(0, 1000)
plt.ylim(1e-5, 1e-2)
plt.legend()
plt.grid()

# save as png image
filename = "EJvsT" + "_mu" + str(mu) + ".png"
plt.savefig(filename)
