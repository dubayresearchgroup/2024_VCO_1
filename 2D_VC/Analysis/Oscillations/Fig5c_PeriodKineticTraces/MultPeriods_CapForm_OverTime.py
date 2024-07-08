import numpy as np
import scipy
import time
import matplotlib.pyplot as plt
from countCapsids import Capsids
import matplotlib.colors



epsA = [1.35, 1.55, 1.75]

period = [4, 2000, 6000, 10000, 15000, 20000]

amp = [0.4]

sigma = 0.25
atoms = 2250
perTri = 15
damp = 0.35

lo = 0.0
hi = 25.48566
partL = 30_000_000
timestep = 0.005
dump= 2000
interval = 50
snaps = int(int(partL / dump) / interval)

head = 9
box1D = 15

dimL = hi - lo
limit = 0.85
tris = atoms / perTri
maxCap = int(tris / 6)
trials = [1, 2, 3]

trialAvg = np.zeros((len(trials),snaps))

PercentCaps = np.zeros((len(period),snaps))

steps = (np.linspace(0, partL, snaps)) 
stepsTau = steps * 0.005
stepsScale = (np.linspace(0, partL, snaps) * 0.005) / 10000

steps_mid = np.linspace(0, partL, snaps) 
t2 = steps_mid * timestep
t2 /= 10000

fig, ax = plt.subplots(3, 1, figsize=(5,15))


colors = ['k', '#3d0f75',  '#fef008', '#70bab5' ,'#3366cc','#0040d9', '#00029e']=

plt.rc('font', size = 22)
plt.rc('axes', labelsize = 28) 

plt.rc('xtick', labelsize = 28) 
plt.rc('ytick', labelsize = 28) 
start = time.time()

trialAvg_stat = np.zeros((len(trials), snaps))


for value in range(len(epsA)):
    for tris in range(len(trials)):

        file = '/Volumes/Niblo/Research/2D_VC/Gaussian_RandomNumber/Equilibrium/0.1/trial%i/zj_real_%.2fpt1.lammpstrj' % (trials[tris],  epsA[value])
        print(file, trials[tris])
        capsids = Capsids(file, head, atoms, snaps, interval, box1D, dimL, limit)
        trialAvg[tris][:] = (100* capsids) / maxCap
       

    PercentCaps_stat = np.mean(trialAvg, axis = 0)
    
    print(PercentCaps_stat)

    ax[value].plot(t2, PercentCaps_stat, zorder =200, color = colors[0], linewidth=2.5)


for a in range(len(amp)):
    for value in range(len(epsA)):
        for p in range(len(period)):
            print(period[p])
            for tris in range(len(trials)):
                eps = epsA[value] 
                file = '/Volumes/Niblo/Research/2D_VC/Gaussian_RandomNumber/Oscillations/Square_Wave/0.1/Amp%.1f/%i/trial%i/zj_real_%.2f_%.1fpt1_%i.lammpstrj' % (amp[a], period[p], trials[tris], eps, amp[a], period[p])
                capsids = Capsids(file, head, atoms, snaps, interval, box1D, dimL, limit)
                trialAvg[tris][:] = (100* capsids) / maxCap
            

            PercentCaps[p][:] = np.mean(trialAvg, axis = 0)
            ax[value].plot(t2, PercentCaps[p], zorder =(100 - (period[p]*0.005)), color = colors[p+1], linewidth=2.5)


        string = '$\epsilon_{\mathrm{avg}}$=%.2f$k_\mathrm{B}T$' %(epsA[value])
        ax[value].set_ylabel('%s ' %(string), fontsize=22)



        ax[value].set_xlim(left = 0)
        ax[value].set_xlim(right = 15)
        ax[value].set_ylim(bottom = 0)
        ax[value].set_ylim(top = 100) 
        ax[value].set_yticks([0, 20, 40, 60, 80, 100])
        ax[value].set_xticks([0, 5, 10, 15])

        ax[value].set_xticklabels([0, 5,10, 15], fontsize=28)
        ax[value].set_yticklabels([0, 20, 40, 60, 80, 100], fontsize=28)

    
Per_tau = [0, 0.02, 10,  30,  50, 75, 100]
print(len(Per_tau))
labels = [1, 10, 20, 30, 40, 50, 60]

cmap = matplotlib.colors.ListedColormap(colors)

norm = matplotlib.colors.Normalize(vmin=0, vmax=60)
sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
halfdist = (labels[2] - labels[1])/2.0
boundaries = np.linspace(labels[0] - halfdist, labels[-1] + halfdist, int(len(labels) + 1.0))


cbar_ax=fig.add_axes([0.99, 0.15, 0.05, 0.7])

cbar = plt.colorbar(sm, spacing='proportional', ticks=labels, orientation='vertical', boundaries=boundaries,label='Period (τ)',cax=cbar_ax)
cbar.set_ticklabels(Per_tau)



fig.supxlabel(r'Time / $10^4 \tau$', fontsize=28, y = 0.04) #, y= -0.0560)
fig.supylabel('Percent Capsid Formation (%)', fontsize=28, x = -0.25)

plt.savefig('PeriodComp_less.png', dpi=300,  bbox_inches = 'tight', pad_inches=0.1)

plt.close()