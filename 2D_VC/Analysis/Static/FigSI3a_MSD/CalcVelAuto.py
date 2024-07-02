import numpy as np
from readTraj_2DVelAuto import ReadSnaps
from neighList_2D import BoxAtoms, GetNeighbors
from operator import itemgetter

def Type3_Data(num_raw, typ, x_raw, y_raw, dimL):

    ##Sort data to get tip
    
    AtomNum = num_raw[np.where(typ == 3)]
    x = x_raw[np.where(typ == 3)] * dimL
    y = y_raw[np.where(typ == 3)] * dimL

    return AtomNum, x, y


def CenterOfMass(tris, atoms, num_raw, typ, x_raw, y_raw, dimL):
    #Calculates the Center of Mass for each triangle

    VerticesAtomNum = np.zeros((tris,3))
    VerticesXCoord = np.zeros((tris,3))
    VerticesYCoord = np.zeros((tris,3))

    ## Find tip of each triangle
    Vert1 = num_raw[np.where(typ == 3)]
    Vert1_x = x_raw[np.where(typ == 3)] * dimL
    Vert1_y = y_raw[np.where(typ == 3)] * dimL

    type3Verts = np.column_stack((Vert1, Vert1_x, Vert1_y))
    type3Verts = sorted(type3Verts, key=itemgetter(0))

    ## Get type b data that can be a potetnial vertice
    possNums = num_raw[np.where(typ == 2)]
    possx = x_raw[np.where(typ == 2)] * dimL
    possy = y_raw[np.where(typ == 2)] * dimL
    potType2Verts = np.column_stack((possNums, possx, possy))
    potType2Verts = sorted(potType2Verts, key=itemgetter(0))

    
    distanceTest = np.zeros(6)
    
    ##Calculate distance between type A vertex and all type B atoms 
    ## to determine which are type B vertices 
    for tri in range(tris):
        shifter = tri * 6
        
        VerticesAtomNum[tri][0] = type3Verts[tri][0]
        VerticesXCoord[tri][0] = type3Verts[tri][1]
        VerticesYCoord[tri][0] = type3Verts[tri][2]
        
        for test in range(6):
            distanceTest[test] = ((abs(type3Verts[tri][1]-potType2Verts[test+shifter][1])**2) + (abs(type3Verts[tri][2]-potType2Verts[test+shifter][2])**2))**1/2

        VerticesAtomNum[tri][1] = potType2Verts[np.argmax(distanceTest)][0]
        VerticesXCoord[tri][1] = potType2Verts[np.argmax(distanceTest)][1]
        VerticesYCoord[tri][1] = potType2Verts[np.argmax(distanceTest)][2]

        Vert3Possibles = set(distanceTest)
        Vert3Possibles.remove(max(Vert3Possibles))
        VerticesAtomNum[tri][2] = possNums[np.argmax(Vert3Possibles)+shifter]
        VerticesXCoord[tri][2] = possx[np.argmax(Vert3Possibles)+shifter]
        VerticesYCoord[tri][2] = possy[np.argmax(Vert3Possibles)+shifter]

    return VerticesAtomNum , VerticesXCoord, VerticesYCoord

    
    

        
def CoM_Velocities(file, hi, tris, head, atoms, snaps, interval, dimL, VerticesAtomNum):

    ##Gets average velocity and position for the COM of each triangle 

    readSnaps = ReadSnaps(file, head, atoms, snaps, interval)
    
    CoM_LastSnap_Pos = np.zeros(2)
    CoM_LastSnap_Vel = np.zeros(2)
    CoM_current_Pos = np.zeros(2)
    CoM_current_Vel = np.zeros(2)
    
    VelCor = np.zeros(snaps)
    XVelSnap = np.zeros(snaps)
    YVelSnap = np.zeros(snaps)
    XPosSnap = np.zeros(snaps)
    YPosSnap = np.zeros(snaps)
    for snap in range(snaps):

        
        num_raw, typ, x_raw, y_raw, x_vel, y_vel = next(readSnaps)
        PerShift = hi / 2

        vert1 = np.zeros(4)
        vert2 = np.zeros(4)
        vert3 = np.zeros(4)

        vert1[0] = x_raw[np.where(num_raw == VerticesAtomNum[0][0])] * dimL #x coord
        vert1[1] = y_raw[np.where(num_raw == VerticesAtomNum[0][0])] * dimL #y coord
        vert1[2] = x_vel[np.where(num_raw == VerticesAtomNum[0][0])]  #vx
        vert1[3] = y_vel[np.where(num_raw == VerticesAtomNum[0][0])]  #vy

        vert2[0] = x_raw[np.where(num_raw == VerticesAtomNum[0][1])] * dimL
        vert2[1] = y_raw[np.where(num_raw == VerticesAtomNum[0][1])] * dimL
        vert2[2] = x_vel[np.where(num_raw == VerticesAtomNum[0][1])]  #vx
        vert2[3] = y_vel[np.where(num_raw == VerticesAtomNum[0][1])]  #vy

        vert3[0] = x_raw[np.where(num_raw == VerticesAtomNum[0][2])] * dimL
        vert3[1] = y_raw[np.where(num_raw == VerticesAtomNum[0][2])] * dimL
        vert3[2] = x_vel[np.where(num_raw == VerticesAtomNum[0][2])]  #vx
        vert3[3] = y_vel[np.where(num_raw == VerticesAtomNum[0][2])]  #vy

        
        CoM_current_Pos[0] = (vert1[0] + vert2[0] + vert3[0]) / 3
        CoM_current_Pos[1] = (vert1[1] + vert2[1] + vert3[1]) / 3

        XPosSnap[snap] = (vert1[0] + vert2[0] + vert3[0]) / 3
        YPosSnap[snap] = (vert1[1] + vert2[1] + vert3[1]) / 3

        XVelSnap[snap] = (vert1[2] + vert2[2] + vert3[2]) / 3
        YVelSnap[snap] = (vert1[3] + vert2[3] + vert3[3]) / 3


        
    return XPosSnap, YPosSnap, XVelSnap, YVelSnap
        

        
        
                
 


        

                
            
        
    

    
