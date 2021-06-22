import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt("output.txt")
N = np.sqrt(data.shape[0])

N = int(N)

x = np.zeros(N)
y = np.zeros(N)
Z1 = np.zeros((N, N))
Z2 = np.zeros((N, N))
for i in range(N) :
    x[i] = data[i, 0]
    y[i] = data[i*N, 2]
    for j in range(N) :
        Z1[i,j] = data[i*N + j, 5] # 3 = phi 4 = Ex, 5 = Ax, 6 = constraint
        Z2[i,j] = data[i*N + j, 6] # 

X,Y = np.meshgrid(x,y)
#plt.contourf(X,Y,Z1)

plt.imshow(Z1, extent=[0, N, 0, N], origin='lower',
           cmap='RdGy', interpolation='nearest') #, vmin=0, vmax=1.0)
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
           cmap='RdGy', interpolation='nearest') #, vmin=0, vmax=1)
plt.colorbar()
ax = plt.gca();
ax.set_xticks(np.arange(0, N, 4));
ax.set_yticks(np.arange(0, N, 4));
plt.xlim((0,N))
plt.ylim((0,N))
ax.grid()
plt.savefig("output_fig2.png")
