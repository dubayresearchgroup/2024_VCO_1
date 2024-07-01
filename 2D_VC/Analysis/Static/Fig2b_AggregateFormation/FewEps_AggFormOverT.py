import numpy as np
import sys
import time
import matplotlib.pyplot as plt
import matplotlib.colors
from countAgg_overT import Agg
#from matplotlib_dashboard import MatplotlibDashboard

epsA = [0.95, 1.35, 1.75]

trials = [1, 2 , 3]

sigma = 0.25
atoms = 2250
perTri = 15
damp = 0.35


lo = 0.0
hi =  25.48566
partL = 30_000_000
timestep = 0.005
dump= 2000
head = 9 #file header for LAMMPS trajectory

parts = range(1)
interval = 50 #skips every 50 snapshots to make analysis quicker
box1D = 15


snaps = int(int(partL / dump) / interval) #how many snapshots within trajectory file
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
avgPercentOverT = np.zeros((4,snaps))
t2 = steps_mid * timestep

##First epsA value
for value in range(1): 
    print(epsA[0])
    start = time.time()
    for tri in range(len(trials)):
        eps = epsA[0] 
        file = '/Volumes/Niblo/Research/2D_VC/Gaussian_RandomNumber/Equilibrium/0.1/trial%i/zj_real_%.2fpt1.lammpstrj' % (trials[tri], eps)
        capSize_overT, sizesOverT = Agg(file, head, atoms, snaps, interval, box1D, dimL, limitcap, perTri)

        avgPerCap_T[tri] = (capSize_overT / (atoms/perTri) ) * 100
        sizesOverT = (sizesOverT / (atoms/perTri) ) * 100
        
        avg[tri] = sizesOverT
        
        print(epsA[value], trials[tri])

    avgPercentOverT_cap0 = np.mean(avgPerCap_T, axis = 0)
    aggData0 = np.mean(avg, axis = 0)

## Second epsA value
for value in range(1): 
    print(epsA[1])
    start = time.time()
    for tri in range(len(trials)):
        eps = epsA[1]
        file = '/Volumes/Niblo/Research/2D_VC/Gaussian_RandomNumber/Equilibrium/0.1/trial%i/zj_real_%.2fpt1.lammpstrj' % (trials[tri], eps)
        capSize_overT, sizesOverT = Agg(file, head, atoms, snaps, interval, box1D, dimL, limitcap, limit, bins, binsC, binW, perTri)

        avgPerCap_T[tri] = (capSize_overT / (atoms/perTri) ) * 100
        sizesOverT = (sizesOverT / (atoms/perTri) ) * 100
        
        
        avg[tri] = sizesOverT
        
    
        print(epsA[value], trials[tri])

    avgPercentOverT_cap1 = np.mean(avgPerCap_T, axis = 0)
    aggData1 = np.mean(avg, axis = 0)


##Third epsA value
for value in range(1): 
    print(epsA[2])
    start = time.time()
    for tri in range(len(trials)):
        eps = epsA[2] 
        file = '/Volumes/Niblo/Research/2D_VC/Gaussian_RandomNumber/Equilibrium/0.1/trial%i/zj_real_%.2fpt1.lammpstrj' % (trials[tri], eps)
        capSize_overT, sizesOverT = Agg(file, head, atoms, snaps, interval, box1D, dimL, limitcap, limit, bins, binsC, binW, perTri)

       

        avgPerCap_T[tri] = (capSize_overT / (atoms/perTri) ) * 100
        sizesOverT = (sizesOverT / (atoms/perTri) ) * 100
        
        
        avg[tri] = sizesOverT

    avgPercentOverT_cap2 = np.mean(avgPerCap_T, axis = 0)
    aggData2 = np.mean(avg, axis = 0)



                
t2 = t2 /10000

### Plotting section of code
plt.rc('font', size = 22)
plt.rc('axes', labelsize = 28)



fig = plt.figure(figsize=(14.4,4.8)) 
gs = fig.add_gridspec(1, 3, wspace = 0, height_ratios=[3])


ax = fig.add_subplot(111)

ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)


ax.set_yticks([])
ax.set_xticks([0])

(ax1, ax2, ax3) = gs.subplots(sharey='row')

##plot first epsA
string = '$\epsilon$ = %.2f $k_\mathrm{B}T$' %(epsA[0])
ax1.text(0.2, 92, string, fontsize = 20)
ax1.plot(t2, aggData0[0], color = '#003366', zorder=0, label = 'Single Particle')
ax1.plot(t2, aggData0[1], color = '#3399b2', zorder=0, label = '2-5 Particles')
ax1.plot(t2, avgPercentOverT_cap0, color = '#66ffff',zorder=10,  label = 'Capsids')
ax1.plot(t2, aggData0[2], color = '#9980ff', zorder=0, label = '6-10 Particles')
ax1.plot(t2, aggData0[3], color = '#bd26ff', zorder=0, label = '11-20 Particles')

ax1.set_ylim(bottom=0)
ax1.set_xlim(left=0)
ax1.set_ylim(top=100)
ax1.set_xlim(right=15)
   
ax1.set_xlim(right=15)
ax1.set_xticks([0, 5, 10])
ax1.set_yticks([0, 20, 40, 60, 80, 100])
ax1.set_ylabel('Percent of Monomers \nin Aggregate Type')


## plot second epsA
string = '$\epsilon$ = %.2f $k_\mathrm{B}T$' %(epsA[1])
ax2.text(0.2, 92, string, fontsize = 20)
ax2.plot(t2, aggData1[0], color = '#003366', zorder=0, label = 'Single Particle')
ax2.plot(t2, aggData1[1], color = '#3399b2', zorder=0, label = '2-5 Particles')
ax2.plot(t2, avgPercentOverT_cap1, color = '#66ffff',zorder=10,  label = 'Capsids')
ax2.plot(t2, aggData1[2], color = '#9980ff', zorder=0, label = '6-10 Particles')
ax2.plot(t2, aggData1[3], color = '#bd26ff', zorder=0, label = '11-20 Particles')

ax2.set_ylim(bottom=0)
ax2.set_xlim(left=0)
ax2.set_ylim(top=100)
ax2.set_xlim(right=15)
   
ax2.set_xlim(right=15)
ax2.set_xticks([0, 5, 10])

## plot third epsA
string = '$\epsilon$ = %.2f $k_\mathrm{B}T$' %(epsA[2])
ax3.text(0.2, 92, string, fontsize = 20)
ax3.plot(t2, aggData2[0], color = '#003366', zorder=0, label = 'Single Particle')
ax3.plot(t2, aggData2[1], color = '#3399b2', zorder=0, label = '2-5 Particles')
ax3.plot(t2, avgPercentOverT_cap2, color = '#66ffff',zorder=10,  label = 'Capsids')
ax3.plot(t2, aggData2[2], color = '#9980ff', zorder=0, label = '6-10 Particles')
ax3.plot(t2, aggData2[3], color = '#bd26ff', zorder=0, label = '11-20 Particles')

ax3.set_ylim(bottom=0)
ax3.set_xlim(left=0)
ax3.set_ylim(top=100)
ax3.set_xlim(right=15)
   
ax3.set_xlim(right=15)
ax3.set_xticks([0, 5, 10, 15])



ax.set_xlabel(r'Time / $10^4 \tau$')
labels2 = ['1', '2-6 ', 'Capsids', '7-10', '11+']
labels = [1, 10, 20, 30, 40] #, 50, 60, 70, 80, 90, 100, 110]
cmap = matplotlib.colors.ListedColormap(['#003366', '#3399b2', '#66ffff', '#9980ff', '#bd26ff'])

        
norm = matplotlib.colors.Normalize(vmin=0, vmax=110)
sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
halfdist = (labels[2] - labels[1])/2.0
boundaries = np.linspace(labels[0] - halfdist, labels[-1] + halfdist, int(len(labels) + 1.0))

l, b, w, h = ax.get_position().bounds
cbar_ax=fig.add_axes([0.245, -0.175, 0.529, 0.07])
cbar = fig.colorbar(sm, cax=cbar_ax, orientation='horizontal',spacing='proportional', ticks=labels, boundaries=boundaries ,label='Aggregate Size (Monomers)') 
cbar.set_ticklabels(labels2)



plt.savefig('test_5trial_rename.png', dpi=300,  bbox_inches = 'tight', pad_inches=0.1)



plt.show()