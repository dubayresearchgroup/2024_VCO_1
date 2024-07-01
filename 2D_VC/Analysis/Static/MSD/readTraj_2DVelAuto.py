import numpy as np

def ReadLines (filename):

    with open(filename, 'r') as file:

        for line in file:

            yield line

def ReadSnap (lines, atoms, head):

    num = np.empty(atoms)
    typ = np.empty(atoms)
    x = np.empty(atoms)
    y = np.empty(atoms)
    vx = np.empty(atoms)
    vy = np.empty(atoms)

    for line in range(head):

        next(lines)

    for atom in range(atoms):

        line = next(lines)

        IDstr, typS, xS, yS, zS, v_x, v_y, v_z = line.split()
        num[atom] = int(IDstr)
        typ[atom] = int(typS)
        x[atom] = float(xS)
        y[atom] = float(yS)
        vx[atom] = float(v_x)
        vy[atom] = float(v_y)
        

    return num, typ, x, y, vx, vy

def ReadSnaps (file, head, atoms, snaps, interval):

    skip = (interval - 1) * (atoms + head)
    lines = ReadLines(file)

    snap = 0
    while snap < snaps:
        
        yield ReadSnap(lines, atoms, head)
        for line in range(skip):
            next(lines)
            
        snap += 1
        
