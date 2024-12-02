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
ham_data = np.loadtxt('data/Ham_Lineout.dat')
mom_dat = np.loadtxt('data/Mom_Lineout.dat')

# Make the plot 
plt.plot(ham_data[1:])
plt.title('lineout of Ham')
plt.xlabel('distance from centre', fontsize=14)
plt.ylabel('Ham', fontsize=14)
plt.savefig('Ham_Lineout.pdf',dpi=256, bbox_inches='tight',pad_inches = 0.1)
plt.close()