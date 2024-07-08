import numpy as np
from readTraj_2D import ReadSnaps
from neighList_2D import BoxAtoms, GetNeighbors
import math

def type2(atomNum, typ, x_raw, y_raw, dimL):
    ## Filter data to get type 2
    num = atomNum[np.where(typ == 2)]
    x = x_raw[np.where(typ == 2)] * dimL
    y = y_raw[np.where(typ == 2)] * dimL

    return x, y, num

def Type3_Data(atomNum, typ, x_raw, y_raw, dimL):
    #Get attractive tip
    num = atomNum[np.where(typ == 3)]
    x = x_raw[np.where(typ == 3)] * dimL
    y = y_raw[np.where(typ == 3)] * dimL

    return x, y, num

def Dist(p1, p2, dimL):
    #Distance calculation 
    dist = abs(p1 - p2)
    if dist > 0.5 * dimL:
        dist = dimL - dist

    return dist

def Capsids(atoms, neighs, x, y, capAtoms, dimL, limitcap): 
    ##Determine number of assembled capsids
    capBoxes = []
    for atom in atoms:

        if atom in capAtoms: capBoxes.append(atom) ; continue

        nearNeighs = []
        for neigh in neighs:
            #Get each neighbor
            if neigh != atom:
                            
                dx = Dist(x[atom], x[neigh], dimL)
                dy = Dist(y[atom], y[neigh], dimL)
                                
                dist = (dx**2 + dy**2)**0.5
                if dist <= limitcap:
                    #if neighbor within limit, could be capsid
                    nearNeighs.append(neigh)
                                
        if len(nearNeighs) == 5:
            ##if there are 5 neighbor, in capsid
            capAtoms.append(atom)
            capAtoms.extend(nearNeighs)
            capBoxes.append(atom)
          
    return capAtoms, capBoxes 


def Agg(file, head, atoms, snaps, interval, box1D, dimL, limitcap, perTri) :
    
    #Construct number of agg, determine what % monomers within aggregate 

    readSnaps = ReadSnaps(file, head, atoms, snaps, interval)
    numAggs = np.empty(snaps)
    
  
    sizesOverT = np.zeros((4,snaps))
    capSize_overT = np.zeros(snaps)
    tots = 0
    for snap in range(snaps):    
        atomNum, typ, x_raw, y_raw = next(readSnaps)
        trianglesX = []
        trianglesY = []
        x, y, num3 = Type3_Data(atomNum, typ, x_raw, y_raw, dimL)
        x2, y2, num2 = type2(atomNum, typ, x_raw, y_raw, dimL)
        
        ##create arrays containing type B vertices
        for i in range(len(num3)):
            newTrix = []
            newTriy = []
            newTrix.append(x[i])
            newTriy.append(y[i])
            atomNum = num3[i]
            triNum = math.ceil(atomNum/perTri)
            for j in range(len(num2)):
                jNum = math.ceil( num2[j] / perTri)
                if jNum == triNum:
                    dx = Dist(x[i], x2[j], dimL)
                    dy = Dist(y[i], y2[j], dimL)
                                    
                    dist = (dx**2 + dy**2)**0.5
            
                    if round(dist,1) == 1.0:
                        newTrix.append(x2[j])
                        newTriy.append(y2[j])
                   

            trianglesX.append(newTrix)
            trianglesY.append(newTriy)

        boxes = BoxAtoms(x, y, box1D, dimL)
        Neighs = GetNeighbors(boxes, box1D)

        number = 0
        
        capAtoms = []
        aggs = []
        
        for box, atoms in boxes.items():
            
            neighs = next(Neighs)
            
            capAtoms, capBoxes = Capsids(atoms, neighs, x, y, capAtoms, dimL, limitcap) ##Get capsids 
            capSize_overT[snap] = len(capAtoms)
            for atom in atoms:

                
                if atom in capBoxes: continue
                if atom in capAtoms: continue

                neighsAtom = []
                
                ##for atoms not in capsids
                for neigh in neighs:
                    if neigh != atom:
                            
                        vertDist = []
                        for i in range(3):
                            for j in range(3):

                                dx = Dist(trianglesX[atom][i], trianglesX[neigh][j], dimL)
                                dy = Dist(trianglesY[atom][i], trianglesY[neigh][j], dimL)
                                vertDist.append((dx**2 + dy**2)**0.5) ##This calculates the distances between all vertices
                        large = max(vertDist)

                        if large <= 2.1:
                            neighsAtom.append(neigh)
                            ##this distance encompasses all distances for different agg. configurations 
                neighsAtom.append(atom)

                exist = 0

                ##combine agg. if they already exist 
                for a in range(len(aggs)):
                    if any(atom in aggs[a] for atom in neighsAtom):
                        for atom in neighsAtom:
                            if atom not in aggs[a]:
                                aggs[a].append(atom)
                        exist = 1
                
                ## if not previously observed,  create new agg 
                if not exist:
                    aggs.append(neighsAtom)

                ##check for duplicates 
                delete = []
                for j in range(len(aggs)):
                    for k in range(len(aggs)):
                        if j < k : 
                            if any(tri in aggs[j] for tri in aggs[k]):
                                for tri in aggs[k]:
                                    if tri not in aggs[j]: 
                                        aggs[j].append(tri)
                                delete.append(k)
                                

                #remove duplicates
                nonDupDelete = []
                for num in delete: 
                    if num not in nonDupDelete: 
                        nonDupDelete.append(num)
       
                nonDupDelete.reverse()                                
                for z in nonDupDelete: del aggs[z]
    
        size_single = []
        size_small = [] #2 - 6 monomers 
        size_medsmall = [] #7 - 10 monomers 
        size_10plus = [] #11 + monomers 

        tots += len(aggs)
        
        ## loop through aggs, place each monomers in array to 
        ## determine % of monomers in agg. type
        for agg in aggs:
            centerAgg = []
            size = len(agg)
           

            for a in range(len(agg)):
                ind = agg[a]
                x_agg = x[ind]
                centerAgg.append(x_agg)
            aggCenter_x = sum(centerAgg) / len(centerAgg)

            if size == 1:
                size_single.append(aggCenter_x)

            elif size >= 2 and size <= 6:
                for a in range(len(agg)):
                    ind = agg[a]
                    x_agg = x[ind]
                    size_small.append(x_agg)
                    

            elif size >= 7 and size <= 10:
                for a in range(len(agg)):
                    ind = agg[a]
                    x_agg = x[ind]
                    size_medsmall.append(x_agg)

            else:
                for a in range(len(agg)):
                    ind = agg[a]
                    x_agg = x[ind]
                    size_10plus.append(x_agg)
                    

        sizeFree = len(size_single)
        size26 = len(size_small)
        size710 = len(size_medsmall)
        size10plus = len(size_10plus)


        sizesOverT[0 , snap] = sizeFree
        sizesOverT[1, snap] = size26
        sizesOverT[2, snap] = size710
        sizesOverT[3, snap] = size10plus
       
    
    
    return capSize_overT, sizesOverT
            
                
       
                
        
    
                    
                        

                                

                            
                            
                            
