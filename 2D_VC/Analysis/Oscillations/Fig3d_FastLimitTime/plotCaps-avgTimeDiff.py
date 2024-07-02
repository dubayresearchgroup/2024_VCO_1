import numpy as np
import scipy
import time
import matplotlib.pyplot as plt
from countCapsids import Capsids
from matplotlib.lines import Line2D
from scipy.interpolate import interp1d
from matplotlib.pyplot import figure
from matplotlib.patches import Patch

start = time.time()
name = 'zj_real'

static_epsA = [0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05, 2.15, 2.25]

epsA = [0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05]
amp = [0.2]
period = [4]

sigma = 0.25
atoms = 2250
perTri = 15
damp = 0.35
repetitions = 1

lo = 0.0
hi = 25.48566
partL = 1500000
timestep = 0.005
dump = 3000
head = 9

numtrials= range(repetitions)
trials = [1,2,3] 
interval = 1
box1D = 15

snaps = int(int(partL / dump) / interval)
dimL = hi - lo
limit = 0.85
tris = atoms / perTri

avgCapsLow = np.zeros((len(static_epsA), len(trials)))
avgCaps = np.zeros((len(static_epsA), len(trials)))
avgCapsHigh = np.zeros((len(static_epsA), len(trials)))

for trial in range(len(trials)):
    for value in range(len(static_epsA)):

        capsidsLow = np.array([])

        file = '/Volumes/Niblo/Research/2D_VC/Gaussian_RandomNumber/Equilibrium/0.1/trial%i/zj_real_%.2fpt1_half.lammpstrj' % (trials[trial],  static_epsA[value])

        partCaps = Capsids(file, head, atoms, snaps, interval, box1D, dimL, limit)
        capsidsLow = np.concatenate((capsidsLow, partCaps), axis = 0)

        avgCapsLow[value][trial] = np.mean(capsidsLow[-snaps:])

        capsids = np.array([])

        file = '/Volumes/Niblo/Research/2D_VC/Gaussian_RandomNumber/Equilibrium/0.1/trial%i/zj_real_%.2fpt1_end.lammpstrj' % (trials[trial],  static_epsA[value])

        partCaps = Capsids(file, head, atoms, snaps, interval, box1D, dimL, limit)
        capsids = np.concatenate((capsids, partCaps), axis = 0)

        avgCaps[value][trial] = np.mean(capsids[-snaps:])

        capsidsHigh = np.array([])

        file = '/Volumes/Niblo/Research/2D_VC/Gaussian_RandomNumber/Equilibrium/0.1/trial%i/zj_real_%.2fpt1_end_long.lammpstrj' % (trials[trial],  static_epsA[value])

        partCaps = Capsids(file, head, atoms, snaps, interval, box1D, dimL, limit)
        capsidsHigh = np.concatenate((capsidsHigh, partCaps), axis = 0)

        avgCapsHigh[value][trial] = np.mean(capsidsHigh[-snaps:])

trialAvg = np.zeros((len(trials), len(epsA)))
maxCap = int(tris / 6)
times = [13500000, 28500000, 118500000]
stdDevO = np.zeros((len(times),len(epsA)))
PercentCaps = np.zeros((len(times),len(epsA)))

for st in range(len(times)):
    for tris in range(len(trials)):
        for value in range(len(epsA)):
            partL = period[0]
            dump = 1
            snaps = int(int(partL / dump) / interval)
           
            capsids = np.array([])

            for seg in range(repetitions):
                print(st, times[st])
                t = times[st] + (30000 * seg)

                file = '/Volumes/Niblo/Research/2D_VC/Gaussian_RandomNumber/Oscillations/Square_Wave/0.1/Amp%.1f/%i/trial%i/zj_real%i_%.2f_%.1fpt1_%i-average%i.lammpstrj' % (amp[0], period[0], trials[tris], seg+1, epsA[value], amp[0], period[0],t)
                print(file, trials[tris])
                
                partCaps = Capsids(file, head, atoms, snaps, interval, box1D, dimL, limit)
                capsids = np.concatenate((capsids, partCaps), axis = 0)
                
            trialAvg[tris][value] = (100* np.mean(capsids)) / maxCap


    PercentCaps[st][:] = np.mean(trialAvg, axis = 0)
    stdDevO[st][:] = np.std(trialAvg, axis = 0)


meanLow = 100 * np.mean(avgCapsLow, axis = 1)  / maxCap
mean = 100 * np.mean(avgCaps, axis = 1)  /maxCap
meanHigh = 100 * np.mean(avgCapsHigh, axis = 1)  / maxCap

stdDevLow = 100 * np.std(avgCapsLow, axis = 1) / maxCap
stdDev = 100 * np.std(avgCaps, axis = 1) / maxCap
stdDevHigh = 100 * np.std(avgCapsHigh, axis = 1) / maxCap



plt.rc('font', size = 18)
plt.rc('axes', labelsize = 22)

handles = [Line2D([0], [0], color='k',label = 'Static'), Line2D([0], [0], color='k', linestyle='dashed', label = '$\pm$0.2 $k_{\mathrm{B}}T$ Osc.'), Patch(facecolor='#6600ff', label=r'75,000 $\tau$'), Patch(facecolor='#990080', label=r'150,000 $\tau$'), Patch(facecolor='#cc0000', label=r'600,000 $\tau$')]


plt.errorbar(static_epsA, meanLow, yerr = stdDevLow, xerr = None, ecolor= '#6600ff', capsize = 5, color = '#6600ff', label = r'75,000 $\tau$')
plt.errorbar(epsA, PercentCaps[0], yerr = stdDevO[0], xerr = None, ecolor= '#6600ff', capsize = 5, color = '#6600ff', linestyle='dashed')

plt.errorbar(static_epsA, mean, yerr = stdDev, xerr = None, ecolor= '#990080', capsize = 5, color = '#990080', label = r'150,000 $\tau$')
plt.errorbar(epsA, PercentCaps[1], yerr = stdDevO[1], xerr = None, ecolor= '#990080', capsize = 5, color = '#990080', linestyle='dashed')

plt.errorbar(static_epsA, meanHigh, yerr = stdDevHigh, xerr = None, ecolor= '#cc0000', capsize = 5, color = '#cc0000', label = r'600,000 $\tau$')
plt.errorbar(epsA, PercentCaps[2], yerr = stdDevO[2], xerr = None, ecolor= '#cc0000', capsize = 5, color = '#cc0000', linestyle='dashed')


plt.yticks([0, 20, 40, 60, 80, 100])
plt.xticks([0.75, 1.0,  1.25, 1.50, 1.75, 2.0, 2.25])
plt.ylim(bottom=0)
plt.ylim(top=100)
plt.xlim(left = 0.75)
plt.xlim(right = 2.25)
plt.legend(handles = handles, prop={'size': 14})#, loc='upper left')

plt.xlabel('Average ε ($k_\mathrm{B}T$)')
plt.savefig('TimeComp-AllTimes_ErrorAllCurves.png', dpi=300, bbox_inches = 'tight', pad_inches=0.1)
plt.close()
  


