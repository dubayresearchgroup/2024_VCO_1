import numpy as np
import matplotlib.pyplot as plt
from CalcVelAuto import CenterOfMass
from CalcVelAuto import CoM_Velocities

import matplotlib.colors
from scipy.optimize import curve_fit

def func(x, n0, t):
    return n0 * np.exp(-x / t)

epsA = [1.15]

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
dump = [10] 
seed = [ 15931, 16409, 29989, 30106, 31400, 40822, 42995, 46622, 50237, 62696, 67443, 68720, 88187, 95059, 98081, 98212, 102290, 879645]


avgSeedsRAC = np.zeros((len(seed),int(segment / dump[0])))


steps = np.linspace(0, segment, int(segment/dump[0]))
time = steps * timestep
zeros = np.zeros(int(segment/dump[0]))

dValues = np.zeros(len(seed))


stagger = 1 #number of snaps to shift
for tri in range(len(seed)):
    resultRAC = np.arange(int(segment / dump[0])).reshape(1,int(segment / dump[0]))

    snaps_seg = int(int(segment / dump[0]) / interval)

    snaps = int(int(partL / dump[0]) / interval)
    file = 'zj_real%.2f_SingleTri_dump%i_seed%i.txt' % (epsA[0], dump[0], seed[tri])
    print(file)

    VerticesAtomNum, VerticesXCoord, VerticesYCoord = CenterOfMass(tris, atoms, num_raw, typ, x_raw, y_raw, dimL)

    XPosSnap, YPosSnap, XVelSnap, YVelSnap, vert1 = CoM_Velocities(file, hi, tris, head, atoms, snaps, interval, dimL, VerticesAtomNum)


    count = 0
    stagStart = 0

    start = stagStart
    stop = snaps_seg

    while start < snaps and stop <= snaps:

        tmp = 0
        seg_RAC = np.zeros(snaps_seg)


        for spot in range(start, stop):
            
            ## Calculate intial vector at t= 0
            xDist = XPosSnap[start] - vert1[start][0] 
            yDist = YPosSnap[start] - vert1[start][1] 
            startV = np.array([xDist, yDist])
            
            mag = np.linalg.norm(startV)
            startV /= mag
            

        

            ##calculate the vector for spot tip to COM at time t
            xDist = XPosSnap[spot] - vert1[spot][0] 
            yDist = YPosSnap[spot] - vert1[spot][1] 
            spotV = np.array([xDist, yDist])
            mag = np.linalg.norm(spotV)
            spotV /= mag

            # #Calculate angle from ref to rotation and current spot 
            seg_RAC[tmp] = np.dot(spotV, startV)  
            tmp += 1

        myRAC = seg_RAC.reshape(1,int(segment / dump[0]))

        resultRAC = np.concatenate((resultRAC, myRAC), axis = 0)

        start += 2
        stop += 2


    resultRAC = np.delete(resultRAC, 0, 0)

    avgSeedsRAC[tri][:] = np.mean(resultRAC, axis = 0)



averageRAC = np.mean(avgSeedsRAC, axis = 0)


tries = 300000
sigma =np.ones(len(averageRAC)) 

popt, pcov = curve_fit(func, time, averageRAC,maxfev = tries, sigma=sigma) 
print(popt, pcov)


fit = np.zeros(len(time))
for x in range(len(time)):
    fit[x] = popt[0] * np.exp(-time[x] / popt[1])
   


plt.rc('font', size = 13)
plt.rc('axes', labelsize = 16)

plt.plot(time, averageRAC, color = 'k',  label = 'Rotational Autocorrelation', zorder=0)
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

