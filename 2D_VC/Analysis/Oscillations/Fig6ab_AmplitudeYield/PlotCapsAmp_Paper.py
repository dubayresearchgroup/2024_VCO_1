import numpy as np
import scipy
import time
import matplotlib.pyplot as plt
from countCapsids import Capsids
import matplotlib.colors


static_epsA = [0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05, 2.15, 2.25, 2.35, 2.45, 2.55, 2.65, 2.75] 


##Amplitudes have different ranges of epsA due to shift and we do not want a negative epsA value
epsA = [1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05, 2.15, 2.25, 2.35, 2.45, 2.55, 2.65, 2.75] 
epsA2 = [0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05, 2.15, 2.25, 2.35, 2.45, 2.55, 2.65, 2.75]


period = [20000]
amp2 = [ 0.2, 0.4]
amp = [ 0.6, 0.8, 1.0]

##To compare against minimum epsA
minepsA2 = [x - amp2[0] for x in epsA2]
minepsA4 = [x - amp2[1] for x in epsA2]
minepsA6 = [x - amp[0] for x in epsA]
minepsA8 = [x - amp[1] for x in epsA]
minepsA10 = [x - amp[2] for x in epsA]


sigma = 0.25
atoms = 2250
perTri = 15
damp = 0.35
repetitions = int((30000000 - 28500000)/30000)
print(repetitions)

lo = 0.0
hi = 25.48566
partL = period[0]
timestep = 0.005

head = 9

interval = 1
box1D = 15

dimL = hi - lo
limit = 0.85
tris = atoms / perTri
maxCap = int(tris / 6)
trials = [1, 2, 3]#, 4, 5]

trialAvg = np.zeros((len(trials), len(epsA)))

stdDev = np.zeros((len(amp),len(epsA)))
PercentCaps = np.zeros((len(amp),len(epsA)))

##Amps 0.6, 0.8, 1.0
for a in range(len(amp)):
    for tris in range(len(trials)):
        for value in range(len(epsA)):
            dump = period[0] / 10
            snaps = int(int(partL / dump) / interval)
            capsids = np.array([])

            for seg in range(repetitions):
                t = 28500000 + (30000 * seg)

                file = '/Volumes/Niblo/Research/2D_VC/Gaussian_RandomNumber/Oscillations/Square_Wave/0.1/Amp%.1f/%i/trial%i/zj_real%i_%.2f_%.1fpt1_%i-average%i.lammpstrj' % (amp[a], period[0], trials[tris], seg+1, epsA[value], amp[a], period[0],t)

                partCaps = Capsids(file, head, atoms, snaps, interval, box1D, dimL, limit)
                capsids = np.concatenate((capsids, partCaps), axis = 0)
                
            trialAvg[tris][value] = (100* np.mean(capsids)) / maxCap
       

    PercentCaps[a][:] = np.mean(trialAvg, axis = 0)
    stdDev[a][:] = np.std(trialAvg, axis = 0)


##Amps 0.2, 0.4
stdDev2 = np.zeros((len(amp2),len(epsA2)))
PercentCaps2 = np.zeros((len(amp2),len(epsA2)))
trialAvg = np.zeros((len(trials), len(epsA2)))
for a in range(len(amp2)):
    for tris in range(len(trials)):
        for value in range(len(epsA2)):
            dump = period[0] / 10
            snaps = int(int(partL / dump) / interval)
            capsids = np.array([])

            for seg in range(repetitions):
                t = 28500000 + (30000 * seg)

                file = '/Volumes/Niblo/Research/2D_VC/Gaussian_RandomNumber/Oscillations/Square_Wave/0.1/Amp%.1f/%i/trial%i/zj_real%i_%.2f_%.1fpt1_%i-average%i.lammpstrj' % (amp2[a], period[0], trials[tris], seg+1, epsA2[value], amp2[a], period[0],t)

                partCaps = Capsids(file, head, atoms, snaps, interval, box1D, dimL, limit)
                capsids = np.concatenate((capsids, partCaps), axis = 0)
                
            trialAvg[tris][value] = (100* np.mean(capsids)) / maxCap
       

    PercentCaps2[a][:] = np.mean(trialAvg, axis = 0)
    stdDev2[a][:] = np.std(trialAvg, axis = 0)



trialAvg_stat = np.zeros((len(trials), len(static_epsA)))

stdDev_stat = np.zeros(len(static_epsA))
PercentCaps_stat = np.zeros(len(static_epsA))
##Static system
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

###PLOTTING 
plt.rc('font', size = 22)
plt.rc('axes', labelsize = 22) 

plt.rc('xtick', labelsize = 18) 
plt.rc('ytick', labelsize = 18) 
fig, ax = plt.subplots(1, 2, figsize=(11, 4.8))


cbar_ax=fig.add_axes([0.93, 0.11, 0.02, 0.77])
amps = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
cmap = matplotlib.colors.ListedColormap(['#660033', '#cc0000', '#de4700', '#e87000', '#f29900', '#ffcc00']) #'#000066', '#2e4c85', '#457394', '#63a6a8', '#82d9bd', '#99ffcc'])
norm = matplotlib.colors.Normalize(vmin=0, vmax=1.0)
sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])

halfdist = (amps[1] - amps[0])/2.0
boundaries = np.linspace(amps[0] - halfdist, amps[-1] + halfdist, int(len(amps) + 1.0))

ax[0].errorbar(static_epsA, PercentCaps_stat, yerr = stdDev_stat, xerr = None, ecolor= '#660033', capsize = 5, color = '#660033')
ax[0].errorbar(epsA2, PercentCaps2[0][:], yerr = None, xerr = None, ecolor= '#cc0000', capsize = 5, color = '#cc0000')
ax[0].errorbar(epsA2, PercentCaps2[1][:], yerr = None, xerr = None, ecolor= '#de4700', capsize = 5, color = '#de4700')
ax[0].errorbar(epsA, PercentCaps[0][:], yerr = None, xerr = None, ecolor= '#e87000', capsize = 5, color = '#e87000')
ax[0].errorbar(epsA, PercentCaps[1][:], yerr = None, xerr = None, ecolor= '#f29900', capsize = 5, color = '#f29900')
ax[0].errorbar(epsA, PercentCaps[2][:], yerr = None, xerr = None, ecolor= '#ffcc00', capsize = 5, color = '#ffcc00')
ax[0].set_xlabel(r'Average $\epsilon$ ($k_{\mathrm{B}}T$)')

ax[0].set_yticks([0, 20, 40, 60, 80, 100])
ax[0].set_ylim(bottom=0)
ax[0].set_ylim(top=100)
ax[0].set_xticks([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
ax[0].set_xlim(left = 0.0)
ax[0].set_xlim(right = 3.0)


ax[1].errorbar(static_epsA, PercentCaps_stat, yerr = stdDev_stat, xerr = None, ecolor= '#660033', capsize = 5, color = '#660033')
ax[1].errorbar(minepsA2, PercentCaps2[0][:], yerr = None, xerr = None, ecolor= '#cc0000', capsize = 5, color = '#cc0000')
ax[1].errorbar(minepsA4, PercentCaps2[1][:], yerr = None, xerr = None, ecolor= '#de4700', capsize = 5, color = '#de4700')
ax[1].errorbar(minepsA6, PercentCaps[0][:], yerr = None, xerr = None, ecolor= '#e87000', capsize = 5, color = '#e87000')
ax[1].errorbar(minepsA8, PercentCaps[1][:], yerr = None, xerr = None, ecolor= '#f29900', capsize = 5, color = '#f29900')
ax[1].errorbar(minepsA10, PercentCaps[2][:], yerr = None, xerr = None, ecolor= '#ffcc00', capsize = 5, color = '#ffcc00')


ax[1].set_yticks([0, 20, 40, 60, 80, 100])
ax[1].set_ylim(bottom=0)
ax[1].set_ylim(top=100)
ax[1].set_xlabel(r'Minimum $\epsilon$ ($k_{\mathrm{B}}T$)')
ax[1].set_xticks([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
ax[1].set_xlim(left = 0.0)
ax[1].set_xlim(right = 3.0)
ax[0].set_ylabel('Percent Capsid Formation')

fig.colorbar(sm, cax=cbar_ax, spacing='proportional', ticks=amps, boundaries=boundaries, format='%2.2g',label=r'Amplitude ($k_{\mathrm{B}}T$)')



plt.savefig('Sq_Per%i_AmpComparison.png' %(period[0]),  dpi = 1000, bbox_inches = 'tight', pad_inches=0.1)
plt.show()






