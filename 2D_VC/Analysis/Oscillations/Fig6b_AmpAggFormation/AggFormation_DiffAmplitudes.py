import numpy as np
import sys
import time
import matplotlib.pyplot as plt
import matplotlib.colors
from countAgg_overT import Agg
# from matplotlib_dashboard import MatplotlibDashboard

epsA = [0.95, 1.25, 1.55]
amp = [0.2, 0.6, 1.0]
period = [20000]
trials = [1, 2, 3]

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
titles = [r'$\epsilon_{\mathrm{min}}$=%.2f$k_\mathrm{B}T$' %(epsA[0]), r'$\epsilon_{\mathrm{min}}$=%.2f$k_\mathrm{B}T$' %(epsA[1]), r'$\epsilon_{\mathrm{min}}$=%.2f$k_\mathrm{B}T$' %(epsA[2])]

snaps = int(int(partL / dump) / interval)

dimL = hi - lo
boxL = hi
limitcap = 0.85
limit = 1.3 
tris = atoms / perTri





avgPerCap_T = np.zeros((len(trials), (snaps * len(parts))))
avg  = np.array([np.zeros((4,(snaps * len(parts)))), np.zeros((4,(snaps * len(parts)))), np.zeros((4,(snaps * len(parts))))]) 

steps_mid = np.linspace(0, partL * len(parts), snaps * len(parts)) 
avgPercentOverT = np.zeros((6,snaps))
t2 = steps_mid * timestep

t2 /= 10000

fig, ax = plt.subplots(3, 3, figsize=(15,15))
plt.rc('font', size = 22)
plt.rc('axes', labelsize = 28) 

plt.rc('xtick', labelsize = 28) 
plt.rc('ytick', labelsize = 28) 
start = time.time()

for value in range(len(epsA)):
    

    for p in range(len(amp)):
        for a in range(len(period)):
            for tri in range(len(trials)):

                eps = epsA[value] + amp[p]
                file = '/Volumes/Niblo/Research/2D_VC/Gaussian_RandomNumber/Oscillations/Square_Wave/0.1/Amp%.1f/%i/trial%i/zj_real_%.2f_%.1fpt1_%i.lammpstrj' % (amp[p], period[a], trials[tri], eps ,amp[p], period[a])
                print(file)
                capSize_overT, sizesOverT = Agg(file, head, atoms, snaps, interval, box1D, dimL, limitcap, limit, perTri)
            

                avgPerCap_T[tri] = (capSize_overT / (atoms/perTri) ) * 100
                sizesOverT = (sizesOverT / (atoms/perTri) ) * 100
                
                
                avg[tri] = sizesOverT
                
            
                print(epsA[value], trials[tri])

    
            avgPercentOverT_cap = np.mean(avgPerCap_T, axis = 0)
            aggData = np.mean(avg, axis = 0)

            if value == 0:
                ax[value, p].set_title(r'Amplitude = %.1f$k_{\mathrm{B}}T$' %(amp[p]), fontsize=22)
            if p == 0: 
                ax[value, p].set_ylabel('%s ' %(titles[value]), fontsize=22)

            string = '$\epsilon_{\mathrm{avg}}$=%.2f$k_\mathrm{B}T$' %(epsA[value] + amp[p])
            ax[value,p].text(0.2, 91, string, zorder = 500, fontsize = 22) #6.25

            string = '$\epsilon_{\mathrm{max}}$=%.2f$k_\mathrm{B}T$' %(epsA[value] + 2 * (amp[p]))
            ax[value,p].text(0.2, 79, string, zorder = 500, fontsize = 22) #6.25



            ax[value, p].set_xlim(left = 0)
            ax[value, p].set_xlim(right = 15)
            ax[value, p].set_ylim(bottom = 0)
            ax[value, p].set_ylim(top = 100) #50.97)
            ax[value, p].set_yticks([0, 20, 40, 60, 80, 100])
            ax[value, p].set_xticks([0, 5, 10, 15])

            ax[value, p].set_xticklabels([0, 5, 10, 15], fontsize=22)
            ax[value, p].set_yticklabels([0, 20, 40, 60, 80, 100], fontsize=22)


            ax[value, p].plot(t2, aggData[0], color = '#003366', label = 'Single Particle')
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

cbar_ax=fig.add_axes([0.93, 0.11, 0.02, 0.77])

cbar = plt.colorbar(sm, cax=cbar_ax, spacing='proportional', orientation='vertical' , ticks=labels, boundaries=boundaries,label='Aggregate Size \n(Monomers)')
cbar.set_ticklabels(labels2)

fig.supxlabel(r'Time / $10^4 \tau$', fontsize=28, y = 0.04) 
fig.supylabel('Monomers in Aggregate Type (%)', fontsize=28, x = 0.025)


plt.savefig('AggComparisonAmplitudes_period%i_cbarVert_v2_rotated_3eps.png' %(period[0]), dpi=300,  bbox_inches = 'tight', pad_inches=0.1)

plt.close()
    
print(time.time() - start)
    
    
