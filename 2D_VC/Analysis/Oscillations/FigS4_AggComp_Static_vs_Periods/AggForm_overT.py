import numpy as np
import sys
import time
import matplotlib.pyplot as plt
import matplotlib.colors
from countAgg_overT import Agg

epsA = [1.15, 1.35, 1.55, 1.75, 1.95]
amp = [0.4]
period = [4, 2000, 20000]
trials = [1, 2, 3 ]

sigma = 0.25
atoms = 2250
perTri = 15
damp = 0.35


lo = 0.0
hi = 25.48566
partL = 30_000_000
timestep = 0.005
dump= 2000
head = 9

parts = range(1)
interval = 50 
box1D = 15
titles = ['Static', r'Period = 0.02$\tau$', r'Period = 10$\tau$',  r'Period = 100$\tau$']


snaps = int(int(partL / dump) / interval)
dimL = hi - lo
boxL = hi
limitcap = 0.85
limit = 1.3 
tris = atoms / perTri


avgCaps = np.zeros(len(epsA))
maxCap = int(tris / 6)

avgPerCap_T = np.zeros((len(trials), (snaps * len(parts))))
avg  = np.array([np.zeros((4,(snaps * len(parts)))), np.zeros((4,(snaps * len(parts)))), np.zeros((4,(snaps * len(parts))))]) 


steps_mid = np.linspace(0, partL * len(parts), snaps * len(parts)) 
avgPercentOverT = np.zeros((6,snaps))
t2 = steps_mid * timestep

t2 /= 10000

fig, ax = plt.subplots(5, 4, figsize=(25,25))
plt.rc('font', size = 24)
plt.rc('axes', labelsize = 26) 

plt.rc('xtick', labelsize = 24) 
plt.rc('ytick', labelsize = 24) 
start = time.time()

for value in range(len(epsA)):
    

    for a in range(len(amp)):
        for p in range(len(period)+1):
            for tri in range(len(trials)):

                if p == 0: 
                    file = '/Volumes/Niblo/Research/2D_VC/Gaussian_RandomNumber/Equilibrium/0.1/trial%i/zj_real_%.2fpt1.lammpstrj' % (trials[tri], epsA[value])
                else:
                    eps = epsA[value] + amp[a]
                    file = '/Volumes/Niblo/Research/2D_VC/Gaussian_RandomNumber/Oscillations/Square_Wave/0.1/Amp%.1f/%i/trial%i/zj_real_%.2f_%.1fpt1_%i.lammpstrj' % (amp[a], period[p-1], trials[tri], epsA[value] ,amp[a], period[p-1])
                print(file)
                capSize_overT, sizesOverT = Agg(file, head, atoms, snaps, interval, box1D, dimL, limitcap, limit, bins, binsC, binW, perTri)
            

                avgPerCap_T[tri] = (capSize_overT / (atoms/perTri) ) * 100
                sizesOverT = (sizesOverT / (atoms/perTri) ) * 100
                
                
                avg[tri] = sizesOverT
                
            
                print(epsA[value], trials[tri])

    
            avgPercentOverT_cap = np.mean(avgPerCap_T, axis = 0)
            aggData = np.mean(avg, axis = 0)

            if value == 0:
                ax[value,p].set_title('%s ' %(titles[p]))

            if p == 0: 
                ax[value,p].set_ylabel(r'$\epsilon_{\mathrm{avg}}$ = %.2f$k_{\mathrm{B}}T$' %(epsA[value]), fontsize=26)

            ax[value,p].set_xlim(left = 0)
            ax[value,p].set_xlim(right = 15)
            ax[value,p].set_ylim(bottom = 0)
            ax[value,p].set_ylim(top = 100) #50.97)
            ax[value,p].set_yticks([0, 20, 40, 60, 80, 100])
            ax[value,p].set_xticks([0, 5, 10, 15])

            ax[value,p].set_xticklabels([0, 5, 10, 15], fontsize=22)
            ax[value,p].set_yticklabels([0, 20, 40, 60, 80, 100], fontsize=22)


            ax[value,p].plot(t2, aggData[0], color = '#003366', label = 'Single Particle')
            ax[value,p].plot(t2, aggData[1], color = '#3399b2', label = '2-5 Particles')
            ax[value,p].plot(t2, aggData[2], color = '#9980ff', label = '6-10 Particles')
            ax[value,p].plot(t2, aggData[3], color = '#bd26ff', label = '11-20 Particles')
            ax[value,p].plot(t2, avgPercentOverT_cap, zorder =200, color = '#66ffff', linewidth=2.5, label = 'Capsids')




                



 
labels2 = ['1', '2-6 ', 'Capsids', '7-10', '11+']
labels = [1, 10, 20, 30, 40]
cmap = matplotlib.colors.ListedColormap(['#003366', '#3399b2', '#66ffff', '#9980ff', '#bd26ff','#ccff33', '#a86da4', '#b8769a', '#c97f8f', '#d98885', '#ea917a', '#fa9a70'])

        
norm = matplotlib.colors.Normalize(vmin=0, vmax=110)
sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
halfdist = (labels[2] - labels[1])/2.0
boundaries = np.linspace(labels[0] - halfdist, labels[-1] + halfdist, int(len(labels) + 1.0))

cbar_ax=fig.add_axes([0.235, 0.045, 0.551, 0.015])
cbar = fig.colorbar(sm, cax=cbar_ax, orientation='horizontal',spacing='proportional', ticks=labels, boundaries=boundaries ,label='Aggregate Size (Monomers)') 
cbar.set_ticklabels(labels2)

fig.supxlabel(r'Time / $10^4 \tau$', y= 0.070)
fig.supylabel('Monomers in Aggregate Type (%)', x = 0.04)


plt.savefig('Agg-DifferentPeriods_v2.png', dpi=600,  bbox_inches = 'tight', pad_inches=0.1)

plt.close()

    
print(time.time() - start)
    
    
