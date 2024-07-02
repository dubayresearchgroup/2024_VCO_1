import numpy as np
import scipy
import time
import matplotlib.pyplot as plt
from countCapsids import Capsids
import matplotlib.colors
from matplotlib import cm


epsA = [0.75, 0.95, 1.15, 1.35, 1.55, 1.75]

period = [20000]
amp = [0.2, 0.4, 0.6,  0.8, 1.0]

sigma = 0.25
atoms = 2250
perTri = 15
damp = 0.35

lo = 0.0
hi = 25.48566
partL = 30_000_000
timestep = 0.005
dump= 2000
interval = 50 
snaps = int(int(partL / dump) / interval)

head = 9
box1D = 15

dimL = hi - lo
limit = 0.85
tris = atoms / perTri
maxCap = int(tris / 6)
trials = [1,2,3] 
trialAvg = np.zeros((len(trials),snaps))


steps = (np.linspace(0, partL, snaps)) 
stepsTau = steps * 0.005
stepsScale = (np.linspace(0, partL, snaps) * 0.005) / 10000

steps_mid = np.linspace(0, partL, snaps) 
t2 = steps_mid * timestep
t2 /= 10000

fig, ax = plt.subplots(2, 3, figsize=(15,10))

plt.rc('font', size = 28)
plt.rc('axes', labelsize = 28) 

plt.rc('xtick', labelsize = 28) 
plt.rc('ytick', labelsize = 28) 
start = time.time()

trialAvg_stat = np.zeros((len(trials), snaps))
epsA_colors = [ '#cc0000', '#de4700', '#e87000', '#f29900', '#ffcc00'] 

for a in range(len(amp)):
    for value in range(len(epsA)):
        for p in range(len(period)):
            print(amp[a], epsA[value])
            for tris in range(len(trials)):
                eps = epsA[value] + amp[a]
                file = '/Volumes/Niblo/Research/2D_VC/Gaussian_RandomNumber/Oscillations/Square_Wave/0.1/Amp%.1f/%i/trial%i/zj_real_%.2f_%.1fpt1_%i.lammpstrj' % (amp[a], period[p], trials[tris], eps, amp[a], period[p])
                capsids = Capsids(file, head, atoms, snaps, interval, box1D, dimL, limit)
                trialAvg[tris][:] = (100* capsids) / maxCap
            

            avgCaps = np.mean(trialAvg, axis = 0)

            row = int(value / 3)
            col = int(value % 3)
            ax[row,col].plot(t2, avgCaps, zorder =200, color = epsA_colors[a], label = '%.1f' %(amp[a])) 
            
            ax[row,col].set_xlim(left = 0)
            ax[row,col].set_xlim(right = 15)
            ax[row,col].set_ylim(bottom = 0)
            ax[row,col].set_ylim(top = 100) #50.97)
            ax[row,col].set_yticks([0, 20, 40, 60, 80, 100])
            ax[row,col].set_xticks([0, 5, 10, 15])

            ax[row,col].set_xticklabels([0, 5, 10, 15], fontsize=22)
            ax[row,col].set_yticklabels([0, 20, 40, 60, 80, 100], fontsize=22)
            string = '$\epsilon_{\mathrm{min}}$ = %.2f $k_\mathrm{B}T$' %(epsA[value])
            ax[row,col].text(0.2, 92, string, fontsize = 20)





cmap = matplotlib.colors.ListedColormap(epsA_colors)
labels = [5, 15, 25, 34, 45]
norm = matplotlib.colors.Normalize(vmin=0, vmax=50)
sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
halfdist = (labels[2] - labels[1])/2.0
boundaries = np.linspace(labels[0] - halfdist, labels[-1] + halfdist, int(len(labels) + 1.0))
cbar_ax=fig.add_axes(([0.93, 0.11, 0.02, 0.77])) #[0.235, 0.025, 0.551, 0.02])



cbar = plt.colorbar(sm, cax=cbar_ax, spacing='proportional', orientation='vertical' , ticks=labels, boundaries=boundaries,label=r'Amplitude ($k_{\mathrm{B}}T$)')
ampLabels = ['0.2', '0.4', '0.6', '0.8', '1.0']
cbar.set_ticklabels(ampLabels)


fig.supxlabel(r'Time / $10^4 \tau$', fontsize=28) #, y= -0.0560)
fig.supylabel('Percent Capsid Formation (%)', fontsize=28, x = 0.05)

plt.savefig('multEps_subplot_three_amps.png', dpi=300,  bbox_inches = 'tight', pad_inches=0.1)

plt.close()