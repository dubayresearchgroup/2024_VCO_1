import numpy as np
import matplotlib.pyplot as plt
from CalcVelAuto import CenterOfMass
from CalcVelAuto import OriginalCoords
from CalcVelAuto import CoM_Velocities

import matplotlib.colors
from scipy.optimize import curve_fit

def func(x, n0, t):
    return n0 * np.exp(-x / t)
    #return a*np.exp(-c*x)+d
    #return m*x+b

o_osc = [0.16, 13.152, 45.344, 72.976, 90.816, 92.992, 92.848, 71.488, 39.376, 17.2, 12.368, 8.4, 3.216]
stdDev = [0.10119289, 1.07164173, 3.22962289, 2.51530992, 2.09373924, 3.49943081, 2.08282885, 3.77159065, 6.40008, 6.17110687, 5.999808, 5.16629074, 3.04101036]
epsA = [1.15] #, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05]

sigma = 0.25
atoms = 15
perTri = 15
damp = 0.35


lo = 0.0
hi = 25.48
partL = 36000
segment = 12000
timestep = 0.005
head = 9

parts = range(1)
interval = 1
box1D = 15


dimL = hi - lo
limit = 0.85
tris = int(atoms / perTri)
trials = [1,2,3,4,5]
dump = [10] #[1, 2, 5, 6, 7, 8, 9, 10,11,12, 13,14, 15,16, 17, 18,19, 20,21, 22, 23, 24, 25,27, 30,31, 33, 35, 38, 40, 41, 43, 48, 45, 50, 55, 60, 65, 70, 75, 80,85,  90, 95, 100, 105, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160,165,166, 170, 175, 180, 185, 190, 195, 200, 205, 215, 225, 230, 235, 240, 245, 250, 255, 260, 265, 275, 280, 285, 290, 295, 300, 312, 325, 330, 335, 340, 345, 350, 355, 360, 365, 370, 375, 380,382, 384, 385, 390, 391, 395, 400, 405, 410, 415, 420, 425, 430, 435, 440, 445, 450, 455, 460, 465, 470, 475, 500]
seed = [ 15931, 16409, 29989, 30106, 31400, 40822, 42995, 46622, 50237, 62696, 67443, 68720, 88187, 95059, 98081, 98212, 102290, 879645]
seed2 = [15931 , 16409, 29989, 31400, 40822, 42995, 46622, 50237, 62696, 67443, 68720, 88187, 95059, 98081, 98212, 879645]

print(len(seed), len(seed2))
#

dumpDiffCoeff = np.zeros(len(dump))
coeff = np.zeros((len(epsA),len(dump)))
combinedVAC =[]
combinedSteps=[]

divider = []

avgSeedsVAC = np.zeros((len(seed),int(segment / dump[0])))
avgSeedsMSD = np.zeros((len(seed),int(segment / dump[0])))

steps = np.linspace(0, segment, int(segment/dump[0]))
time = steps * timestep
zeros = np.zeros(int(segment/dump[0]))

dValues = np.zeros(len(seed))


stagger = 1 #number of snaps to shift
for tri in range(len(seed)):
    resultVAC = np.arange(int(segment / dump[0])).reshape(1,int(segment / dump[0]))
    resultMSD = np.arange(int(segment / dump[0])).reshape(1,int(segment / dump[0]))

    snaps_seg = int(int(segment / dump[0]) / interval)

    snaps = int(int(partL / dump[0]) / interval)
    file = 'zj_real%.2f_SingleTri_dump%i_seed%i.txt' % (epsA[0], dump[0], seed[tri])
    print(file)
    VerticesAtomNum, CoM_initial= OriginalCoords(file, tris, head, atoms, snaps, interval, dimL)

    XPosSnap, YPosSnap, XVelSnap, YVelSnap, vert1 = CoM_Velocities(file, hi, tris, head, atoms, snaps, interval, dimL, VerticesAtomNum)


    count = 0
    stagStart = 0

    start = stagStart
    stop = snaps_seg

    refVector = np.array([1.0,0.0] )
    magnitude_ref = np.linalg.norm(refVector)


    while start < snaps and stop <= snaps:

        tmp = 0
        seg_VAC = np.zeros(snaps_seg)
        seg_MSD = np.zeros(snaps_seg)

        ##calculate the vector for start tip to COM 
        # xDist = XPosSnap[start] - vert1[start][0] #XPosSnap[start]
        # yDist = YPosSnap[start] - vert1[start][1] #YPosSnap[start]
        # startV = np.array([xDist, yDist])
        # magnitude_start = np.linalg.norm(startV)

    
        #calculate angle from ref to rot. at start of window 

        # dot_product = np.dot(startV, refVector)
        # cos_theta_start = dot_product / (magnitude_start * magnitude_ref)
        # angle_radians = np.arccos(np.clip(cos_theta_start, -1.0, 1.0))
        # start_angle_degrees = np.degrees(angle_radians)


        for spot in range(start, stop):

            xDist = XPosSnap[start] - vert1[start][0] #XPosSnap[start]
            yDist = YPosSnap[start] - vert1[start][1] #YPosSnap[start]
            startV = np.array([xDist, yDist])
            
            mag = np.linalg.norm(startV)
            startV /= mag
            

        

            ##calculate the vector for spot tip to COM 
            xDist = XPosSnap[spot] - vert1[spot][0] #XPosSnap[start]
            yDist = YPosSnap[spot] - vert1[spot][1] #YPosSnap[start]
            spotV = np.array([xDist, yDist])
            mag = np.linalg.norm(spotV)
            spotV /= mag


            
        

            dot_product = np.dot(spotV, startV)
          

            # #Calculate angle from ref to rotation and current spot 

            # dot_product = np.dot(spotV, refVector)
            # cos_theta_spot = dot_product / (magnitude_spot * magnitude_ref)
            # angle_radians = np.arccos(np.clip(cos_theta_spot, -1.0, 1.0))
            # spot_angle_degrees = np.degrees(angle_radians)

            seg_VAC[tmp] = np.dot(spotV, startV)  #spot_angle_degrees * start_angle_degrees  #(XVelSnap[spot] * XVelSnap[start]) + ( YVelSnap[spot] * YVelSnap[start])
            #cos_theta_spot * cos_theta_start 
            tmp += 1

        myVAC = seg_VAC.reshape(1,int(segment / dump[0]))

        resultVAC = np.concatenate((resultVAC, myVAC), axis = 0)

        start += 2
        stop += 2


    resultVAC = np.delete(resultVAC, 0, 0)

    print(len(resultVAC))

    avgSeedsVAC[tri][:] = np.mean(resultVAC, axis = 0)

    # tries = 300000
    # sigma =np.ones(len(averagedMSD[20:]))
    # sigma[[0, (len(averagedMSD[20:]))//2, -1]] = 0.01

    # popt, pcov = curve_fit(func, time[20:], averagedMSD[20:], [0.0002,0.0002],maxfev = tries, sigma=sigma)


averageVAC = np.mean(avgSeedsVAC, axis = 0)


tries = 300000
sigma =np.ones(len(averageVAC)) #50:
# sigma[[0, (len(averageMSD[20:]))//2, -1]] = 0.01

popt, pcov = curve_fit(func, time, averageVAC,maxfev = tries, sigma=sigma) 
print(popt, pcov)
# print(popt, epsA[0])
# print('here', time[10], time[20], time[30])

fit = np.zeros(len(time))
for x in range(len(time)):
    fit[x] = popt[0] * np.exp(-time[x] / popt[1])
   


plt.rc('font', size = 13)
plt.rc('axes', labelsize = 16)

plt.plot(time, averageVAC, color = 'k',  label = 'Rotational Autocorrelation', zorder=0)
plt.plot(time, fit, color = '#00ff00', label = 'Exponential Decay', zorder=10)
plt.ylabel(r'$<u(0)u(0+ \tau)$>')
plt.xlabel(r'Time ($\tau$) ')
plt.xlim(left=0)
plt.xlim(right=segment*0.005)
plt.ylim(top=1.2)
plt.ylim(bottom=-0.2)
plt.legend()
plt.axhline(y=0, color='black', linestyle='dashed')
plt.savefig('RotationalAutocorrelation_paper.png', dpi=300, bbox_inches = 'tight', pad_inches=0.1)

plt.show()



# fig, ax = plt.subplots(figsize=[6.4,4.8])
# ax.plot(time, averageMSD, color = 'k',  label = 'MSD', zorder=0)
# ax.plot(time, fit, color = '#00ff00', label = 'Linear Fit', zorder=10)
# ax.axvline(x=0.125, color='blue')
# ax.axvline(x=0.025, color='red')
# ax.set_ylim(bottom=-0.02)
# ax.set_ylim(top=2.0)
# ax.set_xlim(left=0)
# ax.set_xlim(right=20) #20)
# ax.set_ylabel('MSD', color='k')
# ax.set_xlabel(r'Time ($\tau$) ')
# ax.legend(loc='upper left' ) #lower right')

# #axins = ax.inset_axes([0.12, 0.56, 0.4, 0.4])#[0.5, 0.5, 0.47, 0.47])
# axins = ax.inset_axes([0.615, 0.1, 0.35, 0.35])
# axins.plot(time, averageMSD, color='k')
# axins.plot(time, fit, color = '#00ff00', zorder=10)
# x1, x2, y1, y2 = 0, 1, 0, 0.08 #0, 1, 1, 2
# axins.axvline(x=0.025, color='red')
# axins.axvline(x=0.125, color='blue')
# axins.set_xlim(0, 1)
# axins.set_ylim(y1,y2)
# axins.set_xticks([0,  0.5, 1])
# # axins.set_ylabel('MSD', color='k')
# # axins.set_xlabel(r'Time ($\tau$) ')
# axins.set_yticks([0, 0.02, 0.04, 0.06, 0.08])
# plt.savefig('MSD_only_%.2f-fitted_extended.png' %(epsA[0]), dpi=1000, bbox_inches = 'tight', pad_inches=0.1)
