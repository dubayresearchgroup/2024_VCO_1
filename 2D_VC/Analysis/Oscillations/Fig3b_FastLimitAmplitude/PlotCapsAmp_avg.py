import numpy as np
import scipy
import time
import matplotlib.pyplot as plt
from countCapsids import Capsids
import matplotlib.colors


static_epsA = [0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05, 2.15, 2.25]
epsA = [ 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05]


period = [4]
amp = [0.1, 0.2, 0.4]
amp2 = [0.8]

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
trials = [1, 2, 3] 

trialAvg = np.zeros((len(trials), len(epsA)))

stdDev = np.zeros((len(amp),len(epsA)))
PercentCaps = np.zeros((len(amp),len(epsA)))

for a in range(len(amp)):
    for tris in range(len(trials)):
        for value in range(len(epsA)):
            dump = 1
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


plt.rc('font', size = 18)
plt.rc('axes', labelsize = 22)   

plt.errorbar(static_epsA, PercentCaps_stat, yerr = stdDev_stat, xerr = None, ecolor= '#990080', capsize = 5, color = '#990080', label = "Static")
plt.errorbar(epsA, PercentCaps[0][:], yerr = stdDev[0][:], xerr = None, ecolor= '#00cc99', capsize = 5, color = '#00cc99', label = '$\pm$0.1 $k_{\mathrm{B}}T$ Osc.')
plt.errorbar(epsA, PercentCaps[1][:], yerr = stdDev[1][:], xerr = None, ecolor= '#3366b2', capsize = 5, color = '#3366b2', label = '$\pm$0.2 $k_{\mathrm{B}}T$ Osc.')
plt.errorbar(epsA, PercentCaps[2][:], yerr = stdDev[2][:], xerr = None, ecolor= '#6600cc', capsize = 5, color = '#6600cc', label = '$\pm$0.4 $k_{\mathrm{B}}T$ Osc.')


plt.legend(prop={'size': 14})
plt.yticks([0, 20, 40, 60, 80, 100])
plt.ylim(bottom=0)
plt.ylim(top=100)
plt.xlim(right=2.25)
plt.xlim(left=0.75)
plt.xlabel(r'Average $\epsilon$ ($k_{\rm B}T$)')
plt.xticks([0.75, 1.0,  1.25, 1.50, 1.75, 2.0, 2.25])
plt.ylabel('Percent Capsid Formation')
plt.savefig('Sq_Per%i_AmpComparison_Effective_1tris-4steps-ylabel.png' %(period[0]),  dpi = 300, bbox_inches = 'tight', pad_inches=0.1)
plt.show()






å