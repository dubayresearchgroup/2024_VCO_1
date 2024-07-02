import numpy as np
import matplotlib.pyplot as plt
from CalcVelAuto import CenterOfMass
from CalcVelAuto import CoM_Velocities

import matplotlib.colors
from scipy.optimize import curve_fit

def func(x, m, b):
    return m*x+b

epsA = [1.15] 

sigma = 0.25
atoms = 15
perTri = 15
damp = 0.35


lo = 0.0
hi = 25.48
partL = 36000
segment = 4000
timestep = 0.005
head = 9

parts = range(1)
interval = 1
box1D = 15


dimL = hi - lo
limit = 0.85
tris = int(atoms / perTri)

dump = [10] 
seed = [ 15931, 16409, 29989, 30106, 31400, 40822, 42995, 46622, 50237, 62696, 67443, 68720, 88187, 95059, 98081, 98212, 102290, 879645]

avgSeedsVAC = np.zeros((len(seed),int(segment / dump[0])))
avgSeedsMSD = np.zeros((len(seed),int(segment / dump[0])))

steps = np.linspace(0, segment, int(segment/dump[0]))
time = steps * timestep

dValues = np.zeros(len(seed))

stagger = 2 #number of snaps to shift                

for tri in range(len(seed)):

    ##need to initiate with an array of 0s -- this get's deleted later
    resultVAC = np.arange(int(segment / dump[0])).reshape(1,int(segment / dump[0]))
    resultMSD = np.arange(int(segment / dump[0])).reshape(1,int(segment / dump[0]))

    snaps_seg = int(int(segment / dump[0]) / interval)
    
    snaps = int(int(partL / dump[0]) / interval)
    file = 'zj_real%.2f_SingleTri_dump%i_seed%i.txt' % (epsA[0], dump[0], seed[tri])
    print(file)
    

    VerticesAtomNum, VerticesXCoord, VerticesYCoord = CenterOfMass(tris, atoms, num_raw, typ, x_raw, y_raw, dimL)


    XPosSnap, YPosSnap, XVelSnap, YVelSnap = CoM_Velocities(file, hi, tris, head, atoms, snaps, interval, dimL, VerticesAtomNum)

    count = 0
    stagStart = 0

    start = stagStart
    stop = snaps_seg
  
    while start < snaps and stop <= snaps:

        tmp = 0
        ##Storage arrays for each window
        seg_VAC = np.zeros(snaps_seg)
        seg_MSD = np.zeros(snaps_seg)
        
            
        for spot in range(start, stop):

            ## VACF calculation 
            seg_VAC[tmp] = (XVelSnap[spot] * XVelSnap[start]) + ( YVelSnap[spot] * YVelSnap[start])
            
            ## MSD calculation
            xDist = XPosSnap[spot] - XPosSnap[start]
            yDist = YPosSnap[spot] - YPosSnap[start]
            seg_MSD[tmp] = (xDist**2 + yDist**2)

            tmp += 1

        myVAC = seg_VAC.reshape(1,int(segment / dump[0]))
        myMSD = seg_MSD.reshape(1,int(segment / dump[0]))
                        
        resultVAC = np.concatenate((resultVAC, myVAC), axis = 0)
        resultMSD = np.concatenate((resultMSD, myMSD), axis = 0)

        start += 2
        stop += 2

    ##delete intial array of 0s
    resultVAC = np.delete(resultVAC, 0, 0)
    resultMSD = np.delete(resultMSD, 0, 0)

    
    ## Take Average
    averagedMSD = np.mean(resultMSD, axis = 0)
    avgSeedsVAC[tri][:] = np.mean(resultVAC, axis = 0)
    avgSeedsMSD[tri][:] = averagedMSD


    ## Fit data 
    tries = 300000
    sigma =np.ones(len(averagedMSD[20:]))
    sigma[[0, (len(averagedMSD[20:]))//2, -1]] = 0.01

    popt, pcov = curve_fit(func, time[20:], averagedMSD[20:], [0.0002,0.0002],maxfev = tries, sigma=sigma)
    dValues[tri] = (popt[0]/ 4)
    print((popt[0]/ 4))
            
        
## Average data over multiple seeds               
averageVAC = np.mean(avgSeedsVAC, axis = 0)
averageMSD = np.mean(avgSeedsMSD, axis = 0)
errorMSD = np.std(avgSeedsMSD, axis = 0)


print('avg diffusion value', np.mean(dValues), np.std(dValues))



## Fit averaged curve for plotting purposes 
tries = 300000
sigma =np.ones(len(averageMSD[20:]))
sigma[[0, (len(averageMSD[20:]))//2, -1]] = 0.01
                
popt, pcov = curve_fit(func, time[20:], averageMSD[20:], [0.0002,0.0002],maxfev = tries, sigma=sigma)
print(popt, epsA[0])
print('here', time[10], time[20], time[30])    
    
fit = np.zeros(len(time))
for x in range(len(time)):
    fit[x] = (popt[0] * time[x]) + popt[1]

            

plt.rc('font', size = 13)
plt.rc('axes', labelsize = 16)     


fig, ax = plt.subplots(figsize=[6.4,4.8])
ax.plot(time, averageMSD, color = 'k',  label = 'MSD', zorder=0)
ax.plot(time, fit, color = '#00ff00', label = 'Linear Fit', zorder=10)
ax.axvline(x=0.125, color='blue')
ax.axvline(x=0.025, color='red')
ax.set_ylim(bottom=-0.02)
ax.set_ylim(top=2.0)
ax.set_xlim(left=0)
ax.set_xlim(right=20) #20)
ax.set_ylabel('MSD', color='k')
ax.set_xlabel(r'Time ($\tau$) ')
ax.legend(loc='upper left' )

axins = ax.inset_axes([0.615, 0.1, 0.35, 0.35])
axins.plot(time, averageMSD, color='k')
axins.plot(time, fit, color = '#00ff00', zorder=10)
x1, x2, y1, y2 = 0, 1, 0, 0.08 #0, 1, 1, 2
axins.axvline(x=0.025, color='red')
axins.axvline(x=0.125, color='blue')
axins.set_xlim(0, 1)
axins.set_ylim(y1,y2)
axins.set_xticks([0,  0.5, 1])
axins.set_yticks([0, 0.02, 0.04, 0.06, 0.08])

plt.savefig('MSD_only_%.2f-fitted_extended.png' %(epsA[0]), dpi=1000, bbox_inches = 'tight', pad_inches=0.1)

plt.show()


## IF YOU WANT TO PLOT THE VACF: 

# fig, ax1 = plt.subplots()
# ax1.plot(time, averageVAC, '#000099', label = 'C(t)') 
# ax1.plot(time, zeros, 'k')
# ax1.set_xlabel('Simulation Time (τ)')
# ax1.set_ylabel('<V(t)V(0)>', color='#000099')
# ax1.tick_params('y', colors='#000099')
# ax1.set_ylim(bottom=-0.02)
# ax1.set_ylim(top=0.16)
# ax1.set_xlim(left=0)
# ax1.set_xlim(right=20)
# ax1.set_yticks([-0.02, 0.0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16])

# #ax1.set_yticklabels([0,20, 40, 60, 80, 100])

# ax2 = ax1.twinx()
# ax2.plot(time, averageMSD, '#ff0000', label = 'MSD') #, linestyle = '-')
# ax2.plot(0, 0, '#000099', label = 'C(t)') 
# ax2.plot(time[50:], fit[50:], color = '#00ff00', label = 'Linear Fit')
# #ax2.plot(oops, wave, '#f20066', markersize = 0.1, where='post') #,drawstyle = 'steps-mid') #, linestyle = '-')
# ax2.set_ylabel('MSD', color='#ff0000')
# ax2.plot(time, zeros, 'k')
# #ax2.set_ylim([1.2, 2.5])
# ax2.tick_params('y', colors='#ff0000')
# ax2.set_ylim(bottom=-0.25)
# ax2.set_ylim(top=2.25)
# ax2.set_xlim(left=0)
# ax2.set_xlim(right=20)
# ax2.set_yticks([-0.25, 0.0, 0.25, 0.50, 0.75, 1.0, 1.25, 1.50, 1.75, 2.0])
# ax2.legend()



# fig.tight_layout()
# plt.savefig('MSD_VACF_%.2f.png' %(epsA[0]), dpi=300, bbox_inches = 'tight', pad_inches=0.1)
# plt.close()



        
            
                

            



