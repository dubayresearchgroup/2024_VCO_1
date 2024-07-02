import numpy as np
import scipy
import time
import matplotlib.pyplot as plt
from countCapsids import Capsids
import matplotlib.colors
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

static_epsA = [0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05, 2.15, 2.25]

epsA = [0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05]

period = [4]
amp = [0.2]
minepsA = [x - amp[0] for x in epsA]
sigma = 0.25
atoms = 2250
perTri = 15
damp = 0.35
repetitions = int((30000000 - 28500000)/30000)

lo = 0.0
hi = 25.48566


interval = 1

timestep = 0.005

head = 9

interval = 1
box1D = 15

dimL = hi - lo
limit = 0.85
tris = atoms / perTri
maxCap = int(tris / 6)
trials = [1,2,3]

trialAvg = np.zeros((len(trials), len(epsA)))

stdDev_01 = np.zeros((len(period),len(epsA)))
PercentCaps_01 = np.zeros((len(period),len(epsA)))
stdDev_0005 = np.zeros((len(period),len(epsA)))
PercentCaps_0005 = np.zeros((len(period),len(epsA)))

for a in range(len(period)):
    for tris in range(len(trials)):
        for value in range(len(epsA)):
            partL = period[a]
            dump = 1
            snaps = int(int(partL / dump) / interval)
            print(snaps)
            capsids = np.array([])
            dimL = 25.48566
            for seg in range(repetitions):
                t = 28500000 + (30000 * seg)

                file = '/Volumes/Niblo/Research/2D_VC/Gaussian_RandomNumber/Oscillations/Square_Wave/0.1/Amp%.1f/%i/trial%i/zj_real%i_%.2f_%.1fpt1_%i-average%i.lammpstrj' % (amp[0], period[a], trials[tris], seg+1, epsA[value], amp[0], period[a],t)
                print(file, trials[tris])
                
                partCaps = Capsids(file, head, atoms, snaps, interval, box1D, dimL, limit)
                capsids = np.concatenate((capsids, partCaps), axis = 0)
                
            trialAvg[tris][value] = (100* np.mean(capsids)) / maxCap
       

    PercentCaps_01[a][:] = np.mean(trialAvg, axis = 0)
    stdDev_01[a][:] = np.std(trialAvg, axis = 0)

for a in range(len(period)):
    for tris in range(len(trials)):
        for value in range(len(epsA)):
            partL = period[a]
            dump = 1
            snaps = int(int(partL / dump) / interval)
            print(snaps)
            capsids = np.array([])
            
            dimL = 113.975
            for seg in range(repetitions):
                t = 28500000 + (30000 * seg)

                file = '/Volumes/Niblo/Research/2D_VC/Gaussian_RandomNumber/Oscillations/Square_Wave/0.005/Amp%.1f/%i/trial%i/zj_real%i_%.2f_%.1fpt1_%i-average%i.lammpstrj' % (amp[0], period[a], trials[tris], seg+1, epsA[value], amp[0], period[a],t)
                print(file, trials[tris])
                
                partCaps = Capsids(file, head, atoms, snaps, interval, box1D, dimL, limit)
                capsids = np.concatenate((capsids, partCaps), axis = 0)
                
            trialAvg[tris][value] = (100* np.mean(capsids)) / maxCap
       

    PercentCaps_0005[a][:] = np.mean(trialAvg, axis = 0)
    stdDev_0005[a][:] = np.std(trialAvg, axis = 0)
    
trialAvg_stat = np.zeros((len(trials), len(static_epsA)))

stdDev_stat = np.zeros(len(static_epsA))
PercentCaps_stat = np.zeros(len(static_epsA))

for tris in range(len(trials)):
    for value in range(len(static_epsA)):
            dump = 3000
            snaps = int(int(1500000 / dump) / interval)
            capsids = np.array([])

            dimL = 25.48566
    
            file = '/Volumes/Niblo/Research/2D_VC/Gaussian_RandomNumber/Equilibrium/0.1/trial%i/zj_real_%.2fpt1_end.lammpstrj' % (trials[tris],  static_epsA[value])
            print(file, trials[tris])

            partCaps = Capsids(file, head, atoms, snaps, interval, box1D, dimL, limit)
            capsids = np.concatenate((capsids, partCaps), axis = 0)
                
            trialAvg_stat[tris][value] = (100* np.mean(capsids)) / maxCap
       

    PercentCaps_stat[:] = np.mean(trialAvg_stat, axis = 0)
    stdDev_stat[:] = np.std(trialAvg_stat, axis = 0)

print(PercentCaps_stat)
stdDev_stat0005 = np.zeros(len(epsA))
PercentCaps_stat0005 = np.zeros(len(epsA))
trialAvg_stat = np.zeros((len(trials), len(epsA)))
for tris in range(len(trials)):
    for value in range(len(epsA)):
            dump = 3000
            snaps = int(int(1500000 / dump) / interval)
            capsids = np.array([])

            dimL = 113.975
            file = '/Volumes/Niblo/Research/2D_VC/Gaussian_RandomNumber/Equilibrium/0.005/trial%i/zj_real_%.2fpt1_end.lammpstrj' % (trials[tris],  static_epsA[value])
            print(file, trials[tris])

            partCaps = Capsids(file, head, atoms, snaps, interval, box1D, dimL, limit)
            capsids = np.concatenate((capsids, partCaps), axis = 0)
                
            trialAvg_stat[tris][value] = (100* np.mean(capsids)) / maxCap
       

    PercentCaps_stat0005[:] = np.mean(trialAvg_stat, axis = 0)
    stdDev_stat0005[:] = np.std(trialAvg_stat, axis = 0)                      

                    

steps = np.linspace(0, partL, snaps)

t = steps * timestep
        

plt.rc('font', size = 18)
plt.rc('axes', labelsize = 22)


handles = [Line2D([0], [0], color='k',label = 'Static'), Line2D([0], [0], color='k', linestyle='dashed', label = '$\pm$0.2 $k_{\mathrm{B}}T$ Osc.'), Patch(facecolor='#990080', label=r'$\phi$ = 0.1'), Patch(facecolor='#009999', label=r'$\phi$ = 0.005')]

plt.plot([0,0.05], [0,0], color = 'k',label = 'Static' )
plt.plot([0,0.05], [0,0], color = 'k', linestyle='dashed', label = '$\pm$0.2 $k_{B}T$ Osc.' )
plt.errorbar(static_epsA, PercentCaps_stat, yerr = stdDev_stat, xerr = None, ecolor= '#990080', capsize = 5, color = '#990080', label = r'$\phi$ = 0.1')
plt.errorbar(epsA, PercentCaps_01[0][:], yerr = stdDev_01[0], xerr = None, ecolor= '#990080', linestyle='dashed', capsize = 5, color = '#990080')
plt.errorbar(epsA, PercentCaps_stat0005, yerr = stdDev_stat0005, xerr = None, ecolor= '#009999', capsize = 5, color = '#009999', label = r'$\phi$ = 0.005')
plt.errorbar(epsA, PercentCaps_0005[0][:], yerr = stdDev_0005, xerr = None, ecolor= '#009999', linestyle='dashed', capsize = 5, color = '#009999')


plt.legend(handles = handles, prop={'size': 14})
plt.yticks([0, 20, 40, 60, 80, 100])
plt.xticks([0.75, 1.0,  1.25, 1.50, 1.75, 2.0, 2.25])
plt.ylim(top=100)
plt.ylim(bottom=0)
plt.xlim(right = 2.25)
plt.xlim(left = 0.75)

plt.xlabel(r'Average $\epsilon$ ($k_{\rm B}T$)')
plt.ylabel('Percent Capsid Formation')

plt.savefig('ConcComparison_fastlimit.png', dpi=1000, bbox_inches = 'tight', pad_inches=0.1)
plt.close()
