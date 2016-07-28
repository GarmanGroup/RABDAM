

def convertParams(params, margin):
    convParams = [0,0,0,0,0,0]
    convParams[0] = float(params[0]) - margin
    convParams[1] = float(params[1]) + margin
    convParams[2] = float(params[2]) - margin
    convParams[3] = float(params[3]) + margin
    convParams[4] = float(params[4]) - margin
    convParams[5] = float(params[5]) + margin
    return convParams
# end convertParams

def isInXYZparams(atomXYZ, params):
    x = float(atomXYZ[0][0])
    y = float(atomXYZ[1][0])
    z = float(atomXYZ[2][0])
    if float(params[0]) < x < float(params[1]):
        if float(params[2]) < y < float(params[3]):
            if float(params[4]) < z < float(params[5]):
                return True
    return False
# end isInParams
