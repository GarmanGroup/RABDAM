# Copyright Thomas Dixon 2015

def convertParams(params, margin):
    params[0] = float(params[0]) - margin
    params[1] = float(params[1]) + margin
    params[2] = float(params[2]) - margin
    params[3] = float(params[3]) + margin
    params[4] = float(params[4]) - margin
    params[5] = float(params[5]) + margin
    return params    
#end convertParams

def isInXYZparams(atomXYZ, params):
    x = atomXYZ[0][0]
    y = atomXYZ[1][0]
    z = atomXYZ[2][0]
    if params[0] < x < params[1]:
        if params[2] < y < params[3]:
            if params[4] < z < params[5]:
                return True
    return False
#end isInParams