# Copyright Thomas Dixon 2015
def isInXYZparams(atomXYZ, params):
    x = atomXYZ[0]
    y = atomXYZ[1]
    z = atomXYZ[2]
    if atomXYZ[0] < x < atomXYZ[1]:
        if atomXYZ[2] < y < atomXYZ[3]:
            if atomXYZ[4] < z < atomXYZ[5]:
                return True
    return False
#end isInParams