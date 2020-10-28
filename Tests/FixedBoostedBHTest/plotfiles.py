import matplotlib.pyplot as plt
import numpy as np

N = 128
data = np.loadtxt("output.txt")

x = np.zeros(N)
y = np.zeros(N)
Z1 = np.zeros((N, N))
Z2 = np.zeros((N, N))
for i in range(N) :
    x[i] = data[i, 0] - 32
    y[i] = data[i*N, 1] - 32
    for j in range(N) :
        Z1[i,j] = data[i*N + j, 2]
        Z2[i,j] = data[i*N + j, 3]

X,Y = np.meshgrid(x,y)
#plt.contourf(X,Y,Z1)

plt.imshow(Z1, extent=[0, N, 0, N], origin='lower',
           cmap='RdGy', interpolation='nearest')#, vmin=-1e-4, vmax=1e-4)
ax = plt.gca();
ax.set_xticks(np.arange(0, N, 4));
ax.set_yticks(np.arange(0, N, 4));
plt.xlim((0,N))
plt.ylim((0,N))
plt.colorbar()
ax.grid()
plt.savefig("output_fig1.png")

plt.clf()

#plt.contourf(X,Y,Z2)
plt.imshow(Z2, extent=[0, N, 0, N], origin='lower',
           cmap='RdGy', interpolation='nearest')#, vmin=-1e-4, vmax=1e-4)
plt.colorbar()
ax = plt.gca();
ax.set_xticks(np.arange(0, N, 4));
ax.set_yticks(np.arange(0, N, 4));
plt.xlim((0,N))
plt.ylim((0,N))
ax.grid()
plt.savefig("output_fig2.png")
