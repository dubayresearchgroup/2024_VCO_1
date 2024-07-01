import numpy as np
import scipy
import time
import matplotlib.pyplot as plt
from countCapsids import Capsids
from matplotlib.lines import Line2D
from scipy.interpolate import interp1d
from matplotlib.pyplot import figure

start = time.time()
name = 'zj_real'

static_epsA = [0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05, 2.15, 2.25]#, 2.35, 2.45, 2.55]#, 2.65, 2.75]#, 2.85, 2.95, 3.05] #0.75, 0.85, 0.95, 1.05

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

        ##Capsids at 75k tau

        capsidsLow = np.array([])

        file = '/Volumes/Niblo/Research/2D_VC/Gaussian_RandomNumber/Equilibrium/0.1/trial%i/zj_real_%.2fpt1_half.lammpstrj' % (trials[trial],  static_epsA[value])

        partCaps = Capsids(file, head, atoms, snaps, interval, box1D, dimL, limit)
        capsidsLow = np.concatenate((capsidsLow, partCaps), axis = 0)

        avgCapsLow[value][trial] = np.mean(capsidsLow[-snaps:])

        ##Capsids at 150k tau

        capsids = np.array([])

        file = '/Volumes/Niblo/Research/2D_VC/Gaussian_RandomNumber/Equilibrium/0.1/trial%i/zj_real_%.2fpt1_end.lammpstrj' % (trials[trial],  static_epsA[value])

        partCaps = Capsids(file, head, atoms, snaps, interval, box1D, dimL, limit)
        capsids = np.concatenate((capsids, partCaps), axis = 0)

        avgCaps[value][trial] = np.mean(capsids[-snaps:])

        ##Capsids at 600k tau

        capsidsHigh = np.array([])

        file = '/Volumes/Niblo/Research/2D_VC/Gaussian_RandomNumber/Equilibrium/0.1/trial%i/zj_real_%.2fpt1_end_long.lammpstrj' % (trials[trial],  static_epsA[value])

        partCaps = Capsids(file, head, atoms, snaps, interval, box1D, dimL, limit)
        capsidsHigh = np.concatenate((capsidsHigh, partCaps), axis = 0)

        avgCapsHigh[value][trial] = np.mean(capsidsHigh[-snaps:])

##Stats & Std. Deviation        
maxCap = int(tris / 6) 
meanLow = 100 * np.mean(avgCapsLow, axis = 1)  /maxCap
mean = 100 * np.mean(avgCaps, axis = 1)  /maxCap
meanHigh = 100 * np.mean(avgCapsHigh, axis = 1)  /maxCap

stdDevLow = 100 * np.std(avgCapsLow, axis = 1) / maxCap
stdDev = 100 * np.std(avgCaps, axis = 1) / maxCap
stdDevHigh = 100 * np.std(avgCapsHigh, axis = 1) / maxCap

##Legends 
low = Line2D([0], [0], color = '#6600ff', label = r'75,000 $\tau$')
meanL = Line2D([0], [0], color = '#990080', solid_capstyle='round', label = r'150,000 $\tau$')
highL = Line2D([0], [0], color = '#cc0000', label = r'600,000 $\tau$')

##Plotting
plt.rc('font', size = 13)
plt.rc('axes', labelsize = 16)

plt.errorbar(static_epsA, meanLow, yerr = stdDevLow, xerr = None, ecolor= '#6600ff', capsize = 5, color = '#6600ff', label = r'75,000 $\tau$')
plt.errorbar(static_epsA, mean, yerr = stdDev, xerr = None, ecolor= '#990080', capsize = 5, color = '#990080', label = r'150,000 $\tau$')
plt.errorbar(static_epsA, meanHigh, yerr = stdDevHigh, xerr = None, ecolor= '#cc0000', capsize = 5, color = '#cc0000', label = r'600,000 $\tau$')
plt.yticks([0, 20, 40, 60, 80, 100])
plt.xticks([0.75, 1.0,  1.25, 1.50, 1.75, 2.0, 2.25])
plt.ylim(bottom=0)
plt.ylim(top=100)
plt.xlim(left = 0.75)
plt.xlim(right = 2.25)
plt.legend( handles=[low, meanL, highL], loc='upper left')
plt.legend()
plt.xlabel('ε ($k_\mathrm{B}T$)')
plt.ylabel('Percent Capsid Formation')
plt.savefig('TimeComp-AllTimes_trial.png', dpi=300, bbox_inches = 'tight', pad_inches=0.1)
plt.close()
  





