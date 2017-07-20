

def getUnitCellParams(fileName):
    # Parses in the unit cell parameters (the vectors a, b, c, and the angles
    # alpha, beta and gamma) from the input PDB file.

    import math

    fileOpen = open(fileName, 'r')
    for line in fileOpen.readlines():
        if 'CRYST1' in str(line[0:6]):  # Unit cell parameters are stored in this line
            params = line.split()
            a = float(params[1])
            b = float(params[2])
            c = float(params[3])
            alpha = math.radians(float(params[4]))
            beta = math.radians(float(params[5]))
            gamma = math.radians(float(params[6]))
            print 'Unit cell parameters successfully extracted'
            break

    fileOpen.close()
    return (a, b, c, alpha, beta, gamma)


def getAUparams(atomList):
    # Determines the min. and max. values of the x, y and z coordinates in the
    # asymmetric unit. (These are required later when calculating the size
    # of the trimmed atoms box.)

    # Initialises xyz minima and maxima using values from first atom in
    # atomList.
    xMin = float(atomList[0].xyzCoords[0][0])
    xMax = float(atomList[0].xyzCoords[0][0])
    yMin = float(atomList[0].xyzCoords[1][0])
    yMax = float(atomList[0].xyzCoords[1][0])
    zMin = float(atomList[0].xyzCoords[2][0])
    zMax = float(atomList[0].xyzCoords[2][0])

    # Updates xyz minima and maxima as loop through atomList.
    for atm in atomList:
        x = float(atm.xyzCoords[0][0])
        y = float(atm.xyzCoords[1][0])
        z = float(atm.xyzCoords[2][0])

        if x < xMin:
            xMin = x
        elif x > xMax:
            xMax = x
        if y < yMin:
            yMin = y
        elif y > yMax:
            yMax = y
        if z < zMin:
            zMin = z
        elif z > zMax:
            zMax = z

    auParams = [xMin, xMax, yMin, yMax, zMin, zMax]
    return auParams


def convertParams(params, margin):
    # Adds/subtracts the packing density threshold to/from the
    # maximum/minimum x, y and z coordinate values of atoms in the
    # asymmetric unit. The values returned define the boundaries of the
    # trimmed atoms box.

    convParams = [0, 0, 0, 0, 0, 0]
    convParams[0] = float(params[0]) - margin
    convParams[1] = float(params[1]) + margin
    convParams[2] = float(params[2]) - margin
    convParams[3] = float(params[3]) + margin
    convParams[4] = float(params[4]) - margin
    convParams[5] = float(params[5]) + margin
    return convParams


def isInXYZparams(atomXYZ, params):
    # Determines whether the xyz coordinates of atoms in the unit cell 3x3
    # assembly lie within the boundaries of the trimmed atoms box.

    x = float(atomXYZ[0][0])
    y = float(atomXYZ[1][0])
    z = float(atomXYZ[2][0])
    if float(params[0]) < x < float(params[1]):
        if float(params[2]) < y < float(params[3]):
            if float(params[4]) < z < float(params[5]):
                return True
    else:
        return False


def trimAtoms(atomList, params):
    # Removes all atoms with coordinates which lie outside of the trimmed
    # atoms box from the list of atoms in the 3x3 unit cell assembly.
    print 'Excluding the atoms that lie outside of the box'

    trimAtomList = []
    trimAppend = trimAtomList.append

    for atom in atomList:
        atomXYZ = atom.xyzCoords
        if isInXYZparams(atomXYZ, params):
            trimAppend(atom)

    print '%.0f atoms have been retained\n' % int(len(trimAtomList))
    return trimAtomList
