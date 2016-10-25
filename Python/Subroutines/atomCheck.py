

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
