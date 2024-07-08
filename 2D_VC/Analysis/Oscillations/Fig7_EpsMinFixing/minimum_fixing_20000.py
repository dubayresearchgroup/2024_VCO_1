import numpy as np
import scipy
import time
import matplotlib.pyplot as plt
from countCapsids import Capsids
import matplotlib.colors


static_epsA = [0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05, 2.15, 2.25, 2.35, 2.45, 2.55, 2.65, 2.75]#, 2.85, 2.95, 3.05] #0.75, 0.85, 0.95, 1.05

eps_minimum = [1.05, 1.15, 1.25, 1.35, 1.45]


period = [ 20000]

sigma = 0.25
atoms = 2250
perTri = 15
damp = 0.35
repetitions = int((30000000 - 28500000)/30000)

lo = 0.0
hi = 25.48566

timestep = 0.005

head = 9

interval = 1
box1D = 15

dimL = hi - lo
limit = 0.85
tris = atoms / perTri
maxCap = int(tris / 6)
trials = [1,2,3]

plt.rc('font', size = 18)
plt.rc('axes', labelsize = 18) 

plt.rc('xtick', labelsize = 18) 
plt.rc('ytick', labelsize = 18) 
fig, ax = plt.subplots(1, 1)



colors = ['#003366', '#9933ff', '#877aed', '#8099e6', '#73ccd9', '#66ffcc']

## Static Yield Curve ##
trialAvg_stat = np.zeros((len(trials), len(static_epsA)))

stdDev_stat = np.zeros(len(static_epsA))
PercentCaps_stat = np.zeros(len(static_epsA))

for tris in range(len(trials)):
    for value in range(len(static_epsA)):
            dump = 3000
            snaps = int(int(1500000 / dump) / interval)
            capsids = np.array([])

    
            file = '/Volumes/Niblo/Research/2D_VC/Gaussian_RandomNumber/Equilibrium/0.1/trial%i/zj_real_%.2fpt1_end.lammpstrj' % (trials[tris],  static_epsA[value])

            partCaps = Capsids(file, head, atoms, snaps, interval, box1D, dimL, limit)
            capsids = np.concatenate((capsids, partCaps), axis = 0)
                
            trialAvg_stat[tris][value] = (100* np.mean(capsids)) / maxCap
       

    PercentCaps_stat[:] = np.mean(trialAvg_stat, axis = 0)
    stdDev_stat[:] = np.std(trialAvg_stat, axis = 0)

ax.errorbar(static_epsA, PercentCaps_stat, yerr = stdDev_stat, xerr = None, ecolor= colors[0], zorder = 5, capsize = 5, color = colors[0])


##Fixed Minimums
for per in range(len(period)):
    for mini in range(len(eps_minimum)):
        eps_maximum = np.arange(eps_minimum[mini], 3.65, 0.1)

        PercentCaps = np.zeros(len(eps_maximum))
        avg_eps = np.zeros(len(eps_maximum))

        PercentCaps = np.zeros((len(eps_minimum),len(eps_maximum)))
        trialAvg = np.zeros((len(trials), len(eps_maximum)))

        for maxi in range(len(eps_maximum)):
            avg = (eps_minimum[mini] + eps_maximum[maxi]) / 2
            avg_eps[maxi] = avg


            for tris in range(len(trials)):

                partL = period[per]
                dump = period[per] / 10
                snaps = int(int(partL / dump) / interval)
                capsids = np.array([])

                for seg in range(repetitions):
                    t = 28500000 + (30000 * seg)

                    file = '/Volumes/Niblo/Research/2D_VC/Gaussian_RandomNumber/Oscillations/Square_Wave/0.1/Test_Fixing/%i/low_%.2f/trial%i/zj_real_%.2f_min%.2f_max%.2f_%i-average%i.lammpstrj' % (period[per], eps_minimum[mini], trials[tris], avg, eps_minimum[mini], eps_maximum[maxi], period[per],t)
                    print(file)
                    partCaps = Capsids(file, head, atoms, snaps, interval, box1D, dimL, limit)
                    capsids = np.concatenate((capsids, partCaps), axis = 0)

                
                trialAvg[tris][maxi] = (100* np.mean(capsids)) / maxCap

                    
        PercentCaps = np.mean(trialAvg, axis = 0)

        PercentCaps[0] = PercentCaps_stat[3 + mini]

        ax.errorbar(avg_eps, PercentCaps, yerr = None, xerr = None, ecolor= colors[mini+1], zorder = 100, capsize = 5, color = colors[mini+1])
        ax.scatter(avg_eps[0], PercentCaps[0], s=75, color=colors[mini+1], zorder = 500, edgecolors='black')
   




##Finish plotting

cbar_ax=fig.add_axes([0.93, 0.11, 0.02, 0.77])
amps = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
cmap = matplotlib.colors.ListedColormap(colors)



labels = [5, 15, 25, 35, 45, 55]

norm = matplotlib.colors.Normalize(vmin=0, vmax=60)
sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
ax.set_xlabel(r'Average $\epsilon$ ($k_{\mathrm{B}}T$)')
ax.set_xlabel(r'Average $\epsilon$ ($k_{\mathrm{B}}T$)')
ax.set_yticks([0, 20, 40, 60, 80, 100])
ax.set_ylim(bottom=0)
ax.set_ylim(top=100)
ax.set_xticks([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
ax.set_xlim(left = 0.0)
ax.set_xlim(right = 3.0)



ax.set_yticks([0, 20, 40, 60, 80, 100])
ax.set_ylim(bottom=0)
ax.set_ylim(top=100)

ax.set_xticks([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
ax.set_xlim(left = 0.0)
ax.set_xlim(right = 3.0)
ax.set_ylabel('Percent Capsid Formation')

leg = ['Static', '1.05', '1.15', '1.25', '1.35', '1.45']
cbar = fig.colorbar(sm, cax=cbar_ax, spacing='uniform', ticks=labels, label=r'$\epsilon_{\mathrm{min}}$ ($k_{\mathrm{B}}T$)')
cbar.set_ticklabels(leg)


plt.savefig('fixed_minimum_20000k_dot75_include_Start_v2.png',  dpi = 1000, bbox_inches = 'tight', pad_inches=0.1)
plt.show()






