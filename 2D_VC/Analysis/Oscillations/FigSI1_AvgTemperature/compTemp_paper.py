import numpy as np
import matplotlib.pyplot as plt

period = [4, 6, 8, 10, 20, 50, 100, 150, 200,250, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 10000, 15000, 20000, 50000, 100000, 200000, 300000, 500000] 
trials = [1]
epsA =  [1.25]


avgTemp = np.zeros(len(period)+1)

for i in range(len(epsA)):
    avgTemp[0] = epsA[i]
    for p in range(len(period)):

        for t in range(len(trials)):

            file = '../%i/trial%i/TEMP_zj_real%.2f_0.4pt1_%i_end.txt' %(period[p], trials[t], epsA[i], period[p])

            end = np.genfromtxt(file, skip_header=1, usecols=1)

            avgTemp[p+1] = np.round(np.mean(end),6)

            plt.scatter(np.log(period[p]*0.005), np.round(np.mean(end),6), color='k')
            plt.errorbar(np.log(period[p]*0.005), np.round(np.mean(end),6), yerr=np.round(np.std(end),6), color='k')

        print(period[p]*0.005, "&", np.round(np.mean(end),6), "&", np.round(np.std(end),6), "\\\\")



plt.xlim(left=-4)
plt.xlim(right=np.log(5000))
plt.axhline(1.003815, zorder=0, color='lightgrey', lw=90)
plt.axhline(1.003815, zorder=10, color='k',linestyle = '--') # lw=90)
plt.xlabel(r'$\tau_{\rm osc}$') 
plt.xticks(ticks=[-4, -2.3, 0,  2.30258509,  4.60517019, np.log(1000)], labels=[0, 0.10, 1.0,  10,  100, 1000])
plt.ylabel('Observed Temperature')
plt.ylim(bottom=0.90)
plt.ylim(top=1.3)
plt.savefig('temp_v2.png', dpi=600, bbox_inches='tight')
