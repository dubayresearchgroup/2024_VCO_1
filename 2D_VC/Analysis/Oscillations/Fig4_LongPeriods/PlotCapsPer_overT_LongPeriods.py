import numpy as np
import scipy
import time
import matplotlib.pyplot as plt
from countCapsids import Capsids
import matplotlib.colors


epsA = [1.05, 1.35, 1.65] 

period = [10000000] 
amp = [0.4]

sigma = 0.25
atoms = 2250
perTri = 15
damp = 0.35
repetitions = 1

lo = 0.0
hi = 25.48566
partL = 30000000
timestep = 0.005
dump = 2000
interval = 5 
snaps = int(int(partL / dump) / interval)

head = 9
box1D = 15

dimL = hi - lo
limit = 0.85
tris = atoms / perTri
maxCap = int(tris / 6)
trials = [1, 2, 3]

trialAvg = np.zeros((len(trials),snaps))

PercentCaps = np.zeros((len(epsA),snaps))

steps = np.linspace(0, partL, snaps) * 0.005


for a in range(len(period)):
    for value in range(len(epsA)):
        for tris in range(len(trials)):

           
            file = '/Volumes/Niblo/Research/2D_VC/Gaussian_RandomNumber/Oscillations/Square_Wave/0.1/Amp%.1f/%i/trial%i/zj_real_%.2f_%.1fpt1_%i.lammpstrj' % (amp[0], period[a], trials[tris], epsA[value], amp[0], period[a])
            capsids = Capsids(file, head, atoms, snaps, interval, box1D, dimL, limit)
            trialAvg[tris][:] = (100* capsids) / maxCap
        print(tris, period[a])
        
        PercentCaps[value][:] = np.mean(trialAvg, axis = 0)

        

# trials = [1,2, 3]#, 4, 5 ]
# trialAvg = np.zeros((len(trials),snaps))
# static = np.zeros((len(epsA),snaps))
# for value in range(len(epsA)):
#     for tris in range(len(trials)):
#         file = '/Volumes/Niblo/Research/2D_VC/Gaussian_RandomNumber/Equilibrium/0.1/trial%i/zj_real_%.2fpt1.lammpstrj' % (trials[tris], epsA[value])
#         capsids = Capsids(file, head, atoms, snaps, interval, box1D, dimL, limit)
#         trialAvg[tris][:] = (100* capsids) / maxCap

#         print(tris)

#     static[value][:] = np.mean(trialAvg, axis = 0)

low = epsA[0] - amp[0]
high = epsA[0] + amp[0]
steps /= 0.005
rec =[]
for i in range(len(steps)):
    if steps[i] % period[0] < (period[0] / 2):
        #print(steps[i], steps[i] % 20000)
        rec.append(high)
    else:
        rec.append(low)



steps *= 0.005
steps /= 10000 
plt.rc('font', size = 18)
plt.rc('axes', labelsize = 22)       



labels2 = ['0', '6,250', '12,500', '50,000']
labels = [1, 10, 20, 30] 
cmap = matplotlib.colors.ListedColormap(['#000066', '#663da3', '#bf73d9', '#ff99ff', '#ccff33','#ccff33', '#a86da4', '#b8769a', '#c97f8f', '#d98885', '#ea917a', '#fa9a70'])
#cmap = matplotlib.colors.ListedColormap(['#003366', '#734aab', '#ff66ff', '#ffa3c2', '#ffcc99','#ccff99'])

        
norm = matplotlib.colors.Normalize(vmin=0, vmax=110)
sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
halfdist = (labels[2] - labels[1])/2.0
boundaries = np.linspace(labels[0] - halfdist, labels[-1] + halfdist, int(len(labels) + 1.0))


string = r'Period = 50,000 $\tau$'  

fig, ax1 = plt.subplots()
# #cbar = fig.colorbar(sm, spacing='proportional', ticks=labels, boundaries=boundaries,label= r'Period ($\tau$)')
# #cbar.set_ticklabels(labels2)

# #ax1.plot(steps, PercentCaps[0][:], 'k') #, label = '%d', zorder=10) %(period[value])
# #ax1.plot(steps, static, color = '#000066', label = 'Non-Oscillatory')


ax1.set_xlabel(r'Time / $10^4 \tau$')
ax1.text(0.2, 94, string, fontsize = 16)

ax1.plot(steps, PercentCaps[0][:],  color = '#0000cc', zorder=100) 
ax1.plot(steps, PercentCaps[1][:],  color = '#8c3894', zorder=100) 
ax1.plot(steps, PercentCaps[2][:],  color = '#ff6666', zorder=100) 


ax1.set_xticks([0,5,10,15])
ax1.set_ylabel('Percent Capsid Formation', color='k')
ax1.tick_params('y', colors='k')
ax1.set_xlim(left=0)
ax1.set_xlim(right=15)

ax1.set_ylim(bottom=0)
ax1.set_ylim(top=100)


ax1.set_yticklabels([0,20, 40, 60, 80, 100])
ax2 = ax1.twinx()
ax2.plot(steps, rec, 'darkgray', markersize = 0.1, zorder=1, drawstyle = 'steps') #, linestyle = '-')

ax2.set_ylabel(r'$\epsilon$ ($k_{\mathrm{B}}T)$', color='darkgray')
ax2.set_ylim([low-0.065, high+0.065])
s1 = r'$\epsilon_{\mathrm{min}}$'
s2 = r'$\epsilon_{\mathrm{max}}$'
ax2.set_yticks([low, high])
ax2.set_yticklabels([s1, s2])
ax2.tick_params('y', colors='darkgray')

ax1.set_zorder(ax2.get_zorder()+1)
ax1.patch.set_visible(False)
plt.savefig('Formation_overTime_%i_Jan2024-slides1.png' %(period[0]),  dpi = 1000, bbox_inches = 'tight', pad_inches=0.1)
plt.show()




