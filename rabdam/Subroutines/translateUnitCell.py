
# RABDAM
# Copyright (C) 2018 Garman Group, University of Oxford

# This file is part of RABDAM.

# RABDAM is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.

# RABDAM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General
# Public License along with this program.  If not, see
# <http://www.gnu.org/licenses/>.


def convertToCartesian(unit_cell_params):
    # Converts the unit cell parameters into Cartesian vectors.

    import math
    import numpy as np

    print('Converting the unit cell parameters to Cartesian vectors')

    a = unit_cell_params[0]
    b = unit_cell_params[1]
    c = unit_cell_params[2]
    alpha = unit_cell_params[3]
    beta = unit_cell_params[4]
    gamma = unit_cell_params[5]

    # Calculates the volume, v, of a unit parallelepiped with the same angles
    # as the unit cell. The unit cell parameters, plus the volume v, are
    # then used to define the matrix which can convert unit cell into
    # Cartesian vectors.
    v = math.sqrt(1 - math.pow((math.cos(alpha)), 2)
                  - math.pow((math.cos(beta)), 2)
                  - math.pow((math.cos(gamma)), 2)
                  + (2*(math.cos(alpha))*(math.cos(beta))*(math.cos(gamma))))

    a11 = a
    a12 = b*math.cos(gamma)
    a13 = c*math.cos(beta)
    a21 = 0
    a22 = b*math.sin(gamma)
    a23 = c*((math.cos(alpha)-(math.cos(beta)*math.cos(gamma)))/math.sin(gamma))
    a31 = 0
    a32 = 0
    a33 = c*(v/math.sin(gamma))
    conversionMatrix = np.array([[a11, a12, a13],
                                 [a21, a22, a23],
                                 [a31, a32, a33]])

    # The fractional basis vectors of the unit cell are defined, and then are
    # multiplied by the conversion matrix to convert them into their Cartesian
    # counterparts.
    aVector = np.array([[1], [0], [0]])
    bVector = np.array([[0], [1], [0]])
    cVector = np.array([[0], [0], [1]])

    aCartesianVector = np.dot(conversionMatrix, aVector)
    bCartesianVector = np.dot(conversionMatrix, bVector)
    cCartesianVector = np.dot(conversionMatrix, cVector)
    cartesianVectors = (aCartesianVector, bCartesianVector, cCartesianVector)

    print('Conversion complete\n')
    return cartesianVectors


def translateUnitCell(ucAtomList, transAtomList, transAtomIDList,
                      cartesianVectors, aTrans, bTrans, cTrans, count,
                      createAUCpdb, createTApdb):
    # Translates the unit cell +/- 1 units in all dimensions (a, b and c) to
    # generate a 3x3 assembly.

    import numpy as np

    # Converts input unit vectors into matrices.
    aTransMat = np.array([[aTrans], [aTrans], [aTrans]])
    bTransMat = np.array([[bTrans], [bTrans], [bTrans]])
    cTransMat = np.array([[cTrans], [cTrans], [cTrans]])

    # Multiplies the input unit vector matrices by the unit cell Cartesian
    # vectors to determine the distances in x, y and z by which every atom
    # in the unit cell is to be translated. The resultant vectors are then
    # combined into a single translation vector.
    aVec = np.multiply(aTransMat, cartesianVectors[0])
    bVec = np.multiply(bTransMat, cartesianVectors[1])
    cVec = np.multiply(cTransMat, cartesianVectors[2])
    transVector = np.array([aVec, bVec, cVec])
    transVector = np.sum(transVector, axis=0)

    # Each atom in the unit cell is translated as described by the translation
    # matrix, and its new xyz coordinates stored in a list of all atoms in
    # the 3x3 assembly.
    for item in ucAtomList:
        new_atom_xyz = np.add(item.xyzCoords, transVector)
        transAtomList[count][0] = new_atom_xyz[0][0]
        transAtomList[count][1] = new_atom_xyz[1][0]
        transAtomList[count][2] = new_atom_xyz[2][0]
        if createAUCpdb is True or createTApdb is True:
            transAtomIDList[count] = item.atomNum
        count += 1

    print('Successfully translated by (%2sa,%2sb,%2sc) unit cells' % (
        aTrans, bTrans, cTrans
        ))

    return transAtomList, transAtomIDList, count
