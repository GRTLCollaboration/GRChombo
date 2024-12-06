import numpy as np
#from pylab import savefig
import matplotlib.pyplot as plt
import matplotlib.cm as cm

plt.rcParams.update(plt.rcParamsDefault)
plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=15)
plt.rcParams.update({'font.size': 15})
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'figure.figsize'    :  '6, 4.2'})
plt.rcParams.update({'figure.autolayout': True})

# Read data
x_LR = np.loadtxt('data_LR/constraint_lineout.dat')[0][1:]
Ham_LR = np.loadtxt('data_LR/constraint_lineout.dat')[1][1:]
Mom_LR = np.loadtxt('data_LR/constraint_lineout.dat')[2][1:]

x_HR = np.loadtxt('data_HR/constraint_lineout.dat')[0][1:]
Ham_HR = np.loadtxt('data_HR/constraint_lineout.dat')[1][1:]
Mom_HR = np.loadtxt('data_HR/constraint_lineout.dat')[2][1:]

fig = plt.figure(figsize=(8,6))
ax1 = plt.subplot(211)
ax2 = plt.subplot(212)

N_LR = 64
N_HR = 128
n = 2
factor = (N_HR/N_LR)**n

label1 = r"$N^3=128^3$"
label2 = r"$N^3=64^3$"

# Make the plot 
ax1.plot(x_HR, abs(Ham_HR),c="lime", label=label1)
ax1.plot(x_LR, abs(Ham_LR),c="r", label=label2)
ax1.plot(x_LR, abs(Ham_LR/factor),c="k",ls="--", label="2nd order")

ax2.plot(x_HR, Mom_HR,c="lime")
ax2.plot(x_LR, Mom_LR,c="r")
ax2.plot(x_LR, Mom_LR/factor,c="k",ls="--")

ax1.set_yscale('log')
ax2.set_yscale('log')

ax2.set_xlabel('$x$')
ax1.set_ylabel('Ham')
ax2.set_ylabel('Mom')
ax1.set_xticklabels([])

ax1.set_xlim([x_HR[0],x_HR[-1]])
ax2.set_xlim([x_HR[0],x_HR[-1]])

ax1.legend(loc="best")

plt.savefig('constraint_lineout.png',dpi=256, bbox_inches='tight')
plt.close()