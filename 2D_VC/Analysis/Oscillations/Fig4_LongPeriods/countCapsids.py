import numpy as np
from readTraj_2D import ReadSnaps
from neighList_2D import BoxAtoms, GetNeighbors

def Type3_Data(typ, x_raw, y_raw, dimL):
    ## Sort data to pull tips of triangles 
    x = x_raw[np.where(typ == 3)] * dimL
    y = y_raw[np.where(typ == 3)] * dimL

    return x, y

def Dist(p1, p2, dimL):
    #Distance calculation 
    dist = abs(p1 - p2)
    if dist > 0.5 * dimL:
        dist = dimL - dist

    return dist

def Capsids(file, head, atoms, snaps, interval, box1D, dimL, limit):

    readSnaps = ReadSnaps(file, head, atoms, snaps, interval)
    capsids = np.empty(snaps)
    
    for snap in range(snaps):
            
        typ, x_raw, y_raw = next(readSnaps) #get current snapshot
        x, y = Type3_Data(typ, x_raw, y_raw, dimL)

        boxes = BoxAtoms(x, y, box1D, dimL) #filter current positions into grids
        Neighs = GetNeighbors(boxes, box1D) #get nearby neighbors

        capAtoms = []
        for box, atoms in boxes.items(): #iterate through each box of grid
            
            neighs = next(Neighs) #pull neighboring boxes
            for atom in atoms:

                if atom in capAtoms: continue

                nearNeighs = []
                for neigh in neighs: #Determine if neighbor is within limit to be a capsid
                    if neigh != atom:
                            
                        dx = Dist(x[atom], x[neigh], dimL)
                        dy = Dist(y[atom], y[neigh], dimL)
                                
                        dist = (dx**2 + dy**2)**0.5
                        if dist <= limit:
                            nearNeighs.append(neigh)
                                
                if len(nearNeighs) == 5:
                    #if nearby neighbors = 5, this means triangle is in a capsid.
                    capAtoms.append(atom)
                    capAtoms.extend(nearNeighs)

        capsids[snap] = len(capAtoms) / 6 #account for "double" counting as we loop through each potential capsid

    #print(capsids)

    return capsids
