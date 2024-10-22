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

data = np.loadtxt('data/data_out.dat')
rho_dat = np.loadtxt('data/lineout.dat')
t_data = data[:,0]
chi_mean = data[:,3]
rho_mean = data[:,4]
rho_interp = rho_dat[:,1:]

lna = np.log(1/np.sqrt(chi_mean))

plt.plot(t_data,lna,color=cm.Reds(8./10.,1.), marker = ".")
#plt.title('Scale Factor')
#plt.xlim(0,250)
#plt.legend(ncol=2,fontsize=14, bbox_to_anchor=(0., 0.95, 1., 0.102), loc='lower left')
plt.xlabel(r'$t$', fontsize=14)
plt.ylabel(r'ln$(a)$', fontsize=14)
# plt.xscale('log')
plt.savefig('plot_efold.pdf',dpi=256, bbox_inches='tight',pad_inches = 0.1)
plt.close()

plt.plot(lna,rho_mean,color=cm.Reds(8./10.,1.), marker = ".")
#plt.title('Mean of Energy Density')
#plt.xlim(0,250)
#plt.legend(ncol=2,fontsize=14, bbox_to_anchor=(0., 0.95, 1., 0.102), loc='lower left')
plt.xlabel(r'ln$(a)$', fontsize=14)
plt.ylabel(r'$\rho_{\rm mean}$', fontsize=14)
# plt.xscale('log')
plt.savefig('plot_rho_mean.pdf',dpi=256, bbox_inches='tight',pad_inches = 0.1)
plt.close()

L = 1.
N = 128.
dx = L/N
dt_multiplier = 0.5
dt = dx*dt_multiplier # dt = dx * dt_multiplier = 0.5 * 0.25
t = 1/dt
num_t_step = 5 # lineout every time = num_t_step
t_plot = np.arange(0,len(rho_interp),t*num_t_step)

for it_plot, t_plot in enumerate(t_plot):
    plt.plot(rho_interp[int(t_plot)],color=cm.Reds(8./(1+it_plot*10.),1.),
         label='t = '+ str(t_plot*dt), marker = "o")

plt.title('lineout ' + r'$\rho$' + ' along x-axis from origin to box edge')
plt.legend()
#plt.xlim(0,250)
#plt.legend(ncol=2,fontsize=14, bbox_to_anchor=(0., 0.95, 1., 0.102), loc='lower left')
#plt.xlabel(r'ln$(a)$', fontsize=14)
plt.ylabel(r'$\rho$', fontsize=14)
plt.yscale('log')
plt.savefig('plot_lineout.pdf',dpi=256, bbox_inches='tight',pad_inches = 0.1)
plt.close()