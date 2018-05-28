
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


def getAUparams(atomList):
    # Determines the min. and max. values of the x, y and z coordinates in the
    # asymmetric unit. (These are required later when calculating the size
    # of the trimmed atoms box.)

    # Initialises x, y and z minima and maxima using values from first atom in
    # atomList.
    xMin = atomList[0].xyzCoords[0][0]  # xyzCoords is a list of lists, hence
    # the first slice extracts the first list, the second slice extracts the
    # first number in this (single item) list
    xMax = atomList[0].xyzCoords[0][0]
    yMin = atomList[0].xyzCoords[1][0]
    yMax = atomList[0].xyzCoords[1][0]
    zMin = atomList[0].xyzCoords[2][0]
    zMax = atomList[0].xyzCoords[2][0]

    # Updates x, y and z minima and maxima as loop through atomList
    for atm in atomList:
        x = atm.xyzCoords[0][0]
        y = atm.xyzCoords[1][0]
        z = atm.xyzCoords[2][0]

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


def convertParams(params, pdt):
    # Adds/subtracts the packing density threshold to/from the
    # maximum/minimum x, y and z coordinate values of atoms in the
    # asymmetric unit. The values returned define the boundaries of the
    # trimmed atoms box.

    xMin = params[0] - pdt
    xMax = params[1] + pdt
    yMin = params[2] - pdt
    yMax = params[3] + pdt
    zMin = params[4] - pdt
    zMax = params[5] + pdt

    convParams = [xMin, xMax, yMin, yMax, zMin, zMax]
    return convParams


def isInXYZparams(atomXYZ, params):
    # Determines whether the xyz coordinates of atoms in the unit cell 3x3
    # assembly lie within the boundaries of the trimmed atoms box.

    x = atomXYZ[0][0]
    y = atomXYZ[0][1]
    z = atomXYZ[0][2]
    if params[0] < x < params[1]:
        if params[2] < y < params[3]:
            if params[4] < z < params[5]:
                return True
    else:
        return False


def trimAtoms(atomList, params, atom_id_list, createAUCpdb, createTApdb):
    # Removes all atoms with coordinates which lie outside of the trimmed
    # atoms box from the list of atoms in the 3x3 unit cell assembly.

    import numpy as np

    print('Discarding atoms that lie further than 14 Angstroms from the\n'
          'asymmetric unit')

    trimmedAtomList = np.array([[0, 0, 0]])
    trimmedAtomIDList = []
    for num in range(atomList.shape[0]):
        atomXYZ = np.array([[atomList[num][0],
                             atomList[num][1],
                             atomList[num][2]]])
        if isInXYZparams(atomXYZ, params):
            trimmedAtomList = np.append(trimmedAtomList, atomXYZ, axis=0)
            if createAUCpdb is True or createTApdb is True:
                trimmedAtomIDList.append(atom_id_list[num])
    trimmedAtomList = np.delete(trimmedAtomList, 0, axis=0)

    print('--> %s atoms have been retained' % trimmedAtomList.shape[0])
    return trimmedAtomList, trimmedAtomIDList
