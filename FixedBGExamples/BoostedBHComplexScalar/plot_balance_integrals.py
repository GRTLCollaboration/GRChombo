import numpy as np
from pylab import savefig
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

EMS = np.loadtxt('EnMomSourceIntegrals.dat')
F = np.loadtxt('FluxIntegrals.dat')
timedata = EMS[:,0][1:]
dt = timedata[1] - timedata[0]
E = EMS[:,1][1:]*4.
M = EMS[:,2][1:]*4.
S = EMS[:,3][1:]*4.
M0 = M-M[0]
FMi = F[:,2]
FMo = F[:,4]

FMo_dt = np.zeros_like(timedata)
FMi_dt = np.zeros_like(timedata)
S_dt = np.zeros_like(timedata)
for i, t in enumerate(timedata) :
    if (i > 0) :
        FMo_dt[i] += FMo_dt[i-1] + FMo[i] * dt
        FMi_dt[i] += FMi_dt[i-1] + FMi[i] * dt
        S_dt[i] += S_dt[i-1]+ S[i] * dt

plt.plot(timedata,-M0,color=cm.Reds(8./10.,1.),
         label='$\mathcal{Q}_x$ - $\mathcal{Q}_x(t=0)$')
plt.plot(timedata,FMo_dt,color=cm.Blues(7./10.,1),ls='--',
         label=r'$\int \mathcal{F}_{x,{\rm out}}$  dt')
plt.plot(timedata,FMi_dt,color='grey',ls='--',
         label=r'$\int \mathcal{F}_{x, {\rm in}}$  dt')
plt.plot(timedata,S_dt,color=cm.Greens(7./10.,1),ls='-.',
         label='$\int \mathcal{S}_x$  dt')
plt.plot(timedata,FMo_dt - FMi_dt - S_dt,'k:',lw=2.,
         label=r'$\int (\mathcal{F}_{x,{\rm out}} - \mathcal{F}_{x,{\rm in}} - \mathcal{S}_x)$ dt')

plt.xlim(0,2000)
plt.legend(ncol=2,fontsize=14, bbox_to_anchor=(0., 0.95, 1., 0.102), loc='lower left')
plt.xlabel(r'$t/M$', fontsize=14)
plt.savefig('plots/balance_integrals.pdf',dpi=256, bbox_inches='tight',pad_inches = 0.1)
plt.close()
