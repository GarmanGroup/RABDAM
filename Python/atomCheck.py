# Copyright Thomas Dixon 2015
def convertParams(params, margin):
    params[0][0] = params[0][0] - margin
    params[1][0] = params[1][0] + margin
    params[2][0] = params[2][0] - margin
    params[3][0] = params[3][0] + margin
    params[4][0] = params[4][0] - margin
    params[5][0] = params[5][0] + margin
    return params    
#end convertParams

def isInXYZparams(atomXYZ, params):
    x = atomXYZ[0][0]
    y = atomXYZ[1][0]
    z = atomXYZ[2][0]
    if params[0][0] < x < params[1][0]:
        if params[2][0] < y < params[3][0]:
            if params[4][0] < z < params[5][0]:
                return True
    return False
#end isInParams