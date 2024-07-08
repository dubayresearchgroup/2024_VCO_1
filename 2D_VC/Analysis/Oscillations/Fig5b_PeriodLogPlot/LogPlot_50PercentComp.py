import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

plt.rc('font', size = 28)
plt.rc('axes', labelsize = 28)


def lineFromPoints(P, Q):
 
    rise = Q[1] - P[1]
    run = Q[0] - P[0]    
    m = rise / run
    
    b = P[1] - (m * P[0])

    return (50 - b) / m

def FindClosestIndex(capyield): 

    diffArray = capyield - 50
    closestNeg = -100
    closestNegIndex = -1
    closestPos = 100
    closestPosIndex = -1

    ##Determine which value is right above/below 50\% formation
    for index in range(len(diffArray)):
        if diffArray[index] < 0:
            closestNeg = diffArray[index]
            closestNegIndex = index
            continue
        if diffArray[index] > 0:
            closestPos = diffArray[index]
            closestPosIndex = index
            break


    return closestNegIndex, closestPosIndex


def calc(yieldFile, stdFile, tempFile, period, colors, epsVals):

    for p in range(len(period)):
        
        epsA = np.genfromtxt(yieldFile, usecols=0, skip_header=1)
        capyield = np.genfromtxt(yieldFile, usecols=p+1, skip_header=1)
        deviation = np.genfromtxt(stdFile, usecols=p+1, skip_header=1)
        avgTemp = np.genfromtxt(tempFile, usecols=p+1, skip_header=0)


        closestNegIndex, closestPosIndex = FindClosestIndex(capyield)
        point1 = (epsA[closestNegIndex]/ avgTemp[closestNegIndex], capyield[closestNegIndex])
        point2 = (epsA[closestPosIndex]/ avgTemp[closestPosIndex], capyield[closestPosIndex])
        potentialeps = lineFromPoints(point1, point2)

        epsVals.append(potentialeps)

        ### Calculate error up  one direction
        errorUp = capyield + deviation

       
        closestNegIndex, closestPosIndex = FindClosestIndex(errorUp)
        point1 = (epsA[closestNegIndex]/ avgTemp[closestNegIndex], errorUp[closestNegIndex])
        point2 = (epsA[closestPosIndex]/ avgTemp[closestPosIndex], errorUp[closestPosIndex])
        potentialeps_upperError = lineFromPoints(point1, point2)

    
        ##Calculate error down one direction 
        error_down = capyield - deviation

        closestNegIndex, closestPosIndex = FindClosestIndex(error_down)
        point1 = (epsA[closestNegIndex]/ avgTemp[closestNegIndex], error_down[closestNegIndex])
        point2 = (epsA[closestPosIndex]/ avgTemp[closestPosIndex], error_down[closestPosIndex])
        potentialeps_downError = lineFromPoints(point1, point2)

    

        ##Get error ready for plotting
        err = np.zeros((2,1))
        err[0] = abs(potentialeps-potentialeps_downError)
        err[1] = abs(potentialeps_upperError-potentialeps)

        ## plot certain lines in plot 
        ax.axhline(y=1.0770733147726312, zorder = 0, color='grey', lw=12) ##Static yield bar 

        ax.axhline(y=1.077+0.367, zorder = 0, color='lightgrey', lw=12) ##Max yield that we observe
       

        ax.axhline(y=0.5 * (1.077+0.367 + 1.0770733147726312), zorder = 0, color='black', linestyle='--') ##mid point dashed line

        if period[p] == 4: 

            ax.scatter(np.log(period[p]*0.005), potentialeps, s=170, zorder = 10, color=colors[p]) #plots circle
            ax.errorbar(np.log(period[p]*0.005), potentialeps, yerr=err, color= colors[p]) #plots error 
            
        
        else:
            
            ax.scatter(np.log(period[p]*0.005), potentialeps, facecolors= 'none', s=170, edgecolor = colors[p]) #plots open circle
            ax.errorbar(np.log(period[p]*0.005), potentialeps, yerr=err, color= colors[p]) #plots error

 

            ListShownPeriods = [4, 6 , 50, 150, 200, 250, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 10000, 15000, 20000]
            if period[p] in ListShownPeriods: 
                ## check to see if yield curve is shown in Fig 5A 
                ## if so make it a filled in circle
                ax.scatter(np.log(period[p]*0.005), potentialeps, s=170,zorder = 10, facecolors=colors[p], edgecolors = colors[p])


period = [0, 4, 6, 50, 150, 200, 250, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 10000, 15000, 20000]

occurs = np.zeros(len(period))

fig = plt.figure(figsize=(15.0,3.75)) 
ax = fig.add_subplot(111)


##Yield files to read in, colors that correspond to each period 
filesYield = ['yield.txt', 'yield2.txt', 'yield3.txt', 'yield4.txt', 'yield_slowregime.txt'] 
filesStd = ['StdDev.txt', 'StdDev2.txt', 'StdDev3.txt', 'StdDev4.txt', 'StdDev_slowregime.txt'] 
temp = ['Avg_Temps_Amp0.4.txt', 'Avg_Temps_Amp0.4_2.txt', 'Avg_Temps_Amp0.4_3.txt', 'Avg_Temps_Amp0.4_4.txt' , 'AvgTemps_allPers_slow.txt']

pers = [[4, 6, 8, 10, 20, 30, 50, 100, 150, 200, 250], [ 400, 500, 600, 800, 1000, 1400, 1800, 2000], [2500, 3000, 4000, 5000], [6000, 7000, 8000, 10000, 12000, 15000, 20000], [50000, 100000, 300000, 500000]]

colors = [ ['#3d0f75', '#68389e','#9362c8', '#cc99ff','#bf73bf', '#b85c99', '#b04573', '#a82e4d', '#990000','#cc4524','#e06132'], [ '#ff8a48', '#ff9933', '#ffa52d', '#ffb824','#ffc21f','#ffd416', '#ffea0b','#fef008'], [ '#dbee2d',  '#bfe051', '#a9d66d', '#8dc891'] , ['#70bab5', '#479cab','#006699', '#3366cc', '#1a53d2', '#0040d9', '#00029e'], ['#1435a7', '#245bae', '#3380b5', '#42a6bc', '#52ccc3'], ['#1435a7', '#245bae', '#3380b5', '#42a6bc', '#52ccc3']]


period = [4, 6, 8, 10, 20, 30, 50, 100, 150, 200, 250, 400, 500, 600, 800, 1000, 1400, 1800, 2000, 2500, 3000, 4000, 5000, 6000, 7000, 8000, 10000, 12000, 15000, 20000, 100000, 300000, 500000]
tau = np.log(np.array(period) * 0.005)


epsVals = []

for i in range(len(filesYield)):
    calc(filesYield[i], filesStd[i], temp[i], pers[i], colors[i],  epsVals)

epsV = np.array(epsVals)

##Plotting 
ax.set_xlabel(r' $\tau_{\rm osc}$')
ax.set_ylabel('Scaled $\epsilon_{\mathrm{avg}}$ \n at 50% Capsid \n Formation ($k_{\mathrm{B}}T^{\mathrm{(Obs)}}$)')
ax.set_ylim(top=1.50)
ax.set_ylim(bottom=1.0)
ax.set_xlim(left=-4)
ax.set_xlim(right = np.log(5000))
ax.set_xticks([-4, -2.3, 0,  2.30258509,  4.60517019, np.log(1000)]) 
ax.set_xticklabels([0, 0.10, 1.0,  10,  100, 1000])
ax.set_yticks([1.0, 1.25, 1.5])

plt.savefig('50Percent_amp04_LogPlot.png', dpi=300, bbox_inches = 'tight', pad_inches=0.35)
plt.show()