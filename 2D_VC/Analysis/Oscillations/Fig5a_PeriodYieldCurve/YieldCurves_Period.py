from traceback import print_tb
import numpy as np
import scipy
import time
import matplotlib.pyplot as plt
from countCapsids import Capsids
import matplotlib.colors


Per_tau = [0, 0.02, 0.03, 0.25, 0.75, 1.0, 1.25, 2.5, 5, 10, 15, 20, 25, 30, 35, 40, 50, 75, 100]


static_epsA = [0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05]

epsA = [0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05]
period = [250, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 10000, 15000, 20000]

epsA2 = [0.75, 0.85, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05]
period2 = [4, 6 , 50, 150, 200]

amp = [0.4]
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
trials = [1, 2, 3]

trialAvg = np.zeros((len(trials), len(epsA2)))
stdDev = np.zeros((len(period),len(epsA)))
stdDev2 = np.zeros((len(period2),len(epsA2)))
PercentCaps2 = np.zeros((len(period2),len(epsA2)))

avgTemp1 = np.zeros((len(period),len(epsA)))
avgTemp2 = np.zeros((len(period2),len(epsA2)))

for a in range(len(period2)):
    for tris in range(len(trials)):
        for value in range(len(epsA2)):

            tempFile = 'AvgTemps_allPers_2.txt'
            avgTemp2[a, :] = np.genfromtxt(tempFile, usecols=a+1, skip_header=0)


            partL = period2[a]
            if period2[a] < 10: 
                dump = 1
            else:
                dump = period2[a] / 10
            snaps = int(int(partL / dump) / interval)

            print('snaps', snaps)
            capsids = np.array([])

            for seg in range(repetitions):
                t = 28500000 + (30000 * seg)

                file = '/Volumes/Niblo/Research/2D_VC/Gaussian_RandomNumber/Oscillations/Square_Wave/0.1/Amp%.1f/%i/trial%i/zj_real%i_%.2f_%.1fpt1_%i-average%i.lammpstrj' % (amp[0], period2[a], trials[tris], seg+1, epsA2[value], amp[0], period2[a],t)
                print(file, trials[tris])
                
                partCaps = Capsids(file, head, atoms, snaps, interval, box1D, dimL, limit)
                print(partCaps)
                capsids = np.concatenate((capsids, partCaps), axis = 0)
               
                
            trialAvg[tris][value] = (100* np.mean(capsids)) / maxCap
       

    PercentCaps2[a][:] = np.mean(trialAvg, axis = 0)
    stdDev2[a][:] = np.std(trialAvg, axis = 0)


trialAvg = np.zeros((len(trials), len(epsA)))
stdDev = np.zeros((len(period),len(epsA)))
PercentCaps = np.zeros((len(period),len(epsA)))
for a in range(len(period)):
    for tris in range(len(trials)):
        for value in range(len(epsA)):

            tempFile = 'AvgTemps_allPers_1.txt'
            avgTemp1[a, :] = np.genfromtxt(tempFile, usecols=a+1, skip_header=0)


            partL = period[a]
            dump = period[a] / 10
            snaps = int(int(partL / dump) / interval)
            print(snaps)
            capsids = np.array([])

            for seg in range(repetitions):
                t = 28500000 + (30000 * seg)

                file = '/Volumes/Niblo/Research/2D_VC/Gaussian_RandomNumber/Oscillations/Square_Wave/0.1/Amp%.1f/%i/trial%i/zj_real%i_%.2f_%.1fpt1_%i-average%i.lammpstrj' % (amp[0], period[a], trials[tris], seg+1, epsA[value], amp[0], period[a],t)
                print(file, trials[tris])
                
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
            print(file, trials[tris])

            partCaps = Capsids(file, head, atoms, snaps, interval, box1D, dimL, limit)
            capsids = np.concatenate((capsids, partCaps), axis = 0)
                
            trialAvg_stat[tris][value] = (100* np.mean(capsids)) / maxCap
       

    PercentCaps_stat[:] = np.mean(trialAvg_stat, axis = 0)
    stdDev_stat[:] = np.std(trialAvg_stat, axis = 0)


     

                    

steps = np.linspace(0, partL, snaps)

t = steps * timestep
     

plt.rc('font', size = 13)
plt.rc('axes', labelsize = 14)

Per_tau = [0, 0.02, 0.03, 0.25, 0.75, 1.0, 1.25, 2.5, 5, 10, 15, 20, 25, 30, 35, 40, 50, 75, 100]
print(len(Per_tau))
labels = [1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180]
colors = ['k', '#3d0f75', '#6f2390', '#9a34a7','#cc48c2','#ff5c85', '#e06132', '#ff9933','#ffc21f', '#fef008', '#bfe051', '#a9d66d', '#8dc891', '#70bab5', '#479cab', '#006699','#3366cc','#0040d9', '#00029e']


plt.errorbar(static_epsA, PercentCaps_stat, yerr = stdDev_stat, xerr = None, ecolor= colors[0], capsize = 5, color = colors[0])
plt.errorbar(epsA2/avgTemp2[0], PercentCaps2[0][:], yerr = None, xerr = None, ecolor= '#3d0f75', capsize = 5, color = '#3d0f75',  label = r'%.2f $\tau$' %(period2[0]*0.005))
plt.errorbar(epsA2/avgTemp2[1], PercentCaps2[1][:], yerr = None, xerr = None, ecolor= '#68389e', capsize = 5, color = '#68389e', label = r'%.2f $\tau$' %(period2[1]*0.005)) #6f2390')
plt.errorbar(epsA2/avgTemp2[2], PercentCaps2[2][:], yerr = None, xerr = None, ecolor= '#b04573', capsize = 5, color = '#b04573', label = r'%.2f $\tau$' %(period2[2]*0.005))
plt.errorbar(epsA2/avgTemp2[3], PercentCaps2[3][:], yerr = None, xerr = None, ecolor= '#990000', capsize = 5, color = '#990000', label = r'%.2f $\tau$' %(period2[3]*0.005)) #6f2390')
plt.errorbar(epsA2/avgTemp2[4], PercentCaps2[4][:], yerr = None, xerr = None, ecolor= '#cc4524', capsize = 5, color = '#cc4524', label = r'%.2f $\tau$' %(period2[4]*0.005))




plt.errorbar(epsA/avgTemp1[0], PercentCaps[0][:], yerr = None, xerr = None, ecolor= colors[6], capsize = 5, color = colors[6],label = r'%.2f $\tau$' %(period[0]*0.005)) #6f2390'))
plt.errorbar(epsA/avgTemp1[1], PercentCaps[1][:], yerr = None, xerr = None, ecolor= colors[7], capsize = 5, color = colors[7],label = r'%.2f $\tau$' %(period[1]*0.005)) #6f2390'))


plt.errorbar(epsA/avgTemp1[2], PercentCaps[2][:], yerr = None, xerr = None, ecolor= colors[8], capsize = 5, color = colors[8])
plt.errorbar(epsA/avgTemp1[3], PercentCaps[3][:], yerr = None, xerr = None, ecolor= colors[9], capsize = 5, color = colors[9])
plt.errorbar(epsA/avgTemp1[4], PercentCaps[4][:], yerr = None, xerr = None, ecolor= colors[10], capsize = 5, color = colors[10])
plt.errorbar(epsA/avgTemp1[5], PercentCaps[5][:], yerr = None, xerr = None, ecolor= colors[11], capsize = 5, color = colors[11]) #linewidth = 5
plt.errorbar(epsA/avgTemp1[6], PercentCaps[6][:], yerr = None, xerr = None, ecolor= colors[12], capsize = 5, color = colors[12])

plt.errorbar(epsA/avgTemp1[7], PercentCaps[7][:], yerr = None, xerr = None, ecolor= colors[13], capsize = 5, color = colors[13])
plt.errorbar(epsA/avgTemp1[8], PercentCaps[8][:], yerr = None, xerr = None, ecolor= colors[14], capsize = 5, color = colors[14])
plt.errorbar(epsA/avgTemp1[9], PercentCaps[9][:], yerr = None, xerr = None, ecolor= colors[15], capsize = 5, color = colors[15])
plt.errorbar(epsA/avgTemp1[10], PercentCaps[10][:], yerr = None, xerr = None, ecolor= colors[16], capsize = 5, color = colors[16])
plt.errorbar(epsA/avgTemp1[11], PercentCaps[11][:], yerr = None, xerr = None, ecolor= colors[17], capsize = 5, color = colors[17])
plt.errorbar(epsA/avgTemp1[12], PercentCaps[12][:], yerr = None, xerr = None, ecolor= colors[18], capsize = 5, color = colors[18])

plt.hlines(y=50, xmin=1.0, xmax=1.475, linestyle='dashed', color = 'k', zorder=10) #, 1.0, 1.6, linestyle='dashed', color = 'k')

cmap = matplotlib.colors.ListedColormap(['k', '#3d0f75', '#68389e', '#b04573','#990000','#cc4524', '#e06132', '#ff9933','#ffc21f', '#fef008', '#bfe051', '#a9d66d', '#8dc891', '#70bab5', '#479cab', '#006699','#3366cc','#0040d9', '#00029e'])

norm = matplotlib.colors.Normalize(vmin=0, vmax=180)
sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
halfdist = (labels[2] - labels[1])/2.0
boundaries = np.linspace(labels[0] - halfdist, labels[-1] + halfdist, int(len(labels) + 1.0))

colors = ['k', '#3d0f75', '#6f2390', '#9a34a7','#cc48c2','#ff5c85', '#e06132', '#ff9933','#ffc21f', '#fef008', '#bfe051', '#a9d66d', '#8dc891', '#70bab5', '#479cab', '#006699','#3366cc','#0040d9', '#00029e']

cbar = plt.colorbar(sm, spacing='proportional', ticks=labels, boundaries=boundaries,label='Period (τ)', ax=plt.gca())
cbar.set_ticklabels(Per_tau)


plt.yticks([0, 20, 40, 60, 80, 100])
plt.ylim(top=100)
plt.ylim(bottom=0)
plt.xlim(right = 2.05)

plt.xlabel('Scaled Average ε ($k_\mathrm{B}T^{\mathrm{(obs)}}$)')#'Minimum ε ($k_{B}T$)')
plt.ylabel('Percent Capsid Formation')

plt.savefig('Yield_Amp04-avg-all_scaled--Mar24-slide3_14SizeFont.png', dpi=300, bbox_inches = 'tight', pad_inches=0.1)
plt.show()

