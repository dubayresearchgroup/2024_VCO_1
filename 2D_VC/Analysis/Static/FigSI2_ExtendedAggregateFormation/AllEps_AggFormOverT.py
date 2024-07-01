import numpy as np
import sys
import time
import matplotlib.pyplot as plt
import matplotlib.colors
from countAgg_overT import Agg
# from matplotlib_dashboard import MatplotlibDashboard

epsA = [0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55,1.65, 1.75, 1.85, 1.95, 2.05, 2.15]
amp = [0.4]
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
interval = 50 ##skip some snapshots 
box1D = 15


snaps = int(int(partL / dump) / interval)
dimL = hi - lo
boxL = hi
limitcap = 0.85
tris = atoms / perTri


avgCaps = np.zeros(len(epsA))
maxCap = int(tris / 6)


avgPerCap_T = np.zeros((len(trials), (snaps * len(parts))))
avg  = np.array([np.zeros((4,(snaps * len(parts)))), np.zeros((4,(snaps * len(parts)))), np.zeros((4,(snaps * len(parts))))])

steps_mid = np.linspace(0, partL * len(parts), snaps * len(parts)) 
avgPercentOverT = np.zeros((6,snaps))
t2 = steps_mid * timestep

t2 /= 10000

fig, ax = plt.subplots(5, 3, figsize=(15,25))
plt.rc('font', size = 24)
plt.rc('axes', labelsize = 26) 

plt.rc('xtick', labelsize = 24) 
plt.rc('ytick', labelsize = 24) 
start = time.time()

a = 0
p = 0
for value in range(len(epsA)):
    
        for tri in range(len(trials)):

            file = '/Volumes/Niblo/Research/2D_VC/Gaussian_RandomNumber/Equilibrium/0.1/trial%i/zj_real_%.2fpt1.lammpstrj' % (trials[tri], epsA[value])
            
            capSize_overT, sizesOverT = Agg(file, head, atoms, snaps, interval, box1D, dimL, limitcap,  perTri)
        

            avgPerCap_T[tri] = (capSize_overT / (atoms/perTri) ) * 100
            sizesOverT = (sizesOverT / (atoms/perTri) ) * 100
            
            
            avg[tri] = sizesOverT
            
        
            print(epsA[value], trials[tri])


        avgPercentOverT_cap = np.mean(avgPerCap_T, axis = 0)
        aggData = np.mean(avg, axis = 0)

        ##Plotting    
        string = '$\epsilon$ = %.2f $k_\mathrm{B}T$' %(epsA[value])
        if a > 3:
            ax[a,p].text(0.7, 70, string, zorder = 500, fontsize = 20)
        else:
            ax[a,p].text(0.2, 90, string, zorder = 500, fontsize = 20)


        ax[a,p].set_xlim(left = 0)
        ax[a,p].set_xlim(right = 15)
        ax[a,p].set_ylim(bottom = 0)
        ax[a,p].set_ylim(top = 100) #50.97)
        ax[a,p].set_yticks([0, 20, 40, 60, 80, 100])
        ax[a,p].set_xticks([0, 5, 10, 15])

        ax[a,p].set_xticklabels([0, 5, 10, 15], fontsize=22)
        ax[a,p].set_yticklabels([0, 20, 40, 60, 80, 100], fontsize=22)


        ax[a,p].plot(t2, aggData[0], color = '#003366', label = 'Single Particle')
        ax[a,p].plot(t2, aggData[1], color = '#3399b2', label = '2-5 Particles')
        ax[a,p].plot(t2, aggData[2], color = '#9980ff', label = '6-10 Particles')
        ax[a,p].plot(t2, aggData[3], color = '#bd26ff', label = '11-20 Particles')
        ax[a,p].plot(t2, avgPercentOverT_cap, zorder =200, color = '#66ffff', linewidth=2.5, label = 'Capsids')

        p += 1

        ##Determine which subplot is next
        if p % 3 == 0: 
            p = 0
            a += 1


labels2 = ['1', '2-6 ', 'Capsids', '7-10', '11+']
labels = [1, 10, 20, 30, 40] #, 50, 60, 70, 80, 90, 100, 110]
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


plt.savefig('StaticAggs_v3-Font_smaller.png', dpi=600,  bbox_inches = 'tight', pad_inches=0.1)

plt.close()

    
print(time.time() - start)
    
    
