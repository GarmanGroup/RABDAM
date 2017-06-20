

class atom(object):
    # A class to assign properties to atoms (specifically all properties
    # provided in the input PDB file, plus values relating to the B_damage
    # calculation).

    def __init__(self, lineidentifier='', atomnum=0, residuenum=0,
                 atomtype='', resitype='', chainID='', xyz_coords=[],
                 atomidentifier='', bfactor=0, occupancy=1, charge='',
                 packingdensity=0, avrg_bfactor=0, bdamage=0):
        self.lineID = lineidentifier
        self.atomNum = atomnum
        self.resiNum = residuenum
        self.atomType = atomtype
        self.resiType = resitype
        self.chainID = chainID
        self.xyzCoords = xyz_coords
        self.atomID = atomidentifier
        self.bFactor = bfactor
        self.occupancy = occupancy
        self.charge = charge
        self.pd = packingdensity
        self.avrg_bf = avrg_bfactor
        self.bd = bdamage


def downloadPDB(PDBcode, PDBdirectory, pathToPDB):
    # Downloads and saves PDB file from the RSCB PDB website.

    import os
    import urllib2

    urlText = 'http://www.rcsb.org/pdb/files/%s.pdb' % PDBcode
    os.makedirs(PDBdirectory)
    print 'Directory %s created' % PDBdirectory
    origPDB = urllib2.urlopen(urlText)
    print 'Downloaded PDB file from %s' % urlText
    localFile = open(pathToPDB, 'w')
    localFile.write(origPDB.read())
    print 'PDB file saved to %s' % pathToPDB
    localFile.close()


def copyPDB(pathToPDB, disk, newPathToPDB, PDBdirectory):
    # Copies .pdb / .txt file from any provided file path to Logfiles
    # directory.

    import os

    owd = os.getcwd()
    os.chdir('/')
    os.chdir(disk)
    origPDB = open(pathToPDB, 'r')
    os.chdir(owd)
    os.makedirs(PDBdirectory)
    localFile = open(newPathToPDB, 'w')
    localFile.write(origPDB.read())
    print 'PDB file copied to %s' % newPathToPDB
    localFile.close()
    origPDB.close()


def full_atom_list(fileName):
    # Parses in input PDB file (= PDBCUR output) and returns a list of all
    # atoms in the file with their associated attributes as assigned in the
    # 'atom' class. The lines above and below the 'ATOM'/'HETATM' information
    # in the input PDB file are also stored in lists for use in writing new
    # PDB files.

    # Initialise lists for storage of subsections of input PDB file
    bof = []
    bAppend = bof.append
    atomList = []
    alAppend = atomList.append
    eof = []
    eAppend = eof.append

    beforeAtoms = True  # Lines before the 'ATOM'/'HETATM' information are stored in 'bof' list

    fileOpen = open(fileName, 'r')
    for line in fileOpen.readlines():
        if str((line[0:6]).strip()) in ['ATOM', 'HETATM']:
            beforeAtoms = False  # Lines after the 'ATOM' / 'HETATM' information are stored in 'eof' list
            y = atom()
            y.lineID = str(line[0:6].strip())
            y.atomNum = int(line[6:11].strip())
            y.atomType = str(line[12:16].strip())
            y.resiType = str(line[17:20].strip())
            y.chainID = str(line[21:22].strip())
            y.resiNum = int(line[22:26].strip())
            y.xyzCoords = [[float(line[30:38].strip())],
                           [float(line[38:46].strip())],
                           [float(line[46:54].strip())]]
            y.occupancy = float(line[54:60].strip())
            y.bFactor = float(line[60:66].strip())
            y.atomID = str(line[76:78].strip())
            y.charge = str(line[78:80].strip())
            alAppend(y)

        else:
            if beforeAtoms is True:
                bAppend(line)
            else:
                eAppend(line)

    fileOpen.close()
    print 'Finished reading in atoms --> %d atoms found' % int(len(atomList))
    print 'in %s' % fileName

    return bof, atomList, eof


def b_damage_atom_list(fileName, atomList, HETATM, protOrNA, addAtoms,
                       removeAtoms):
    # Trims a copy of the atomList returned by the 'full_atom_list' function
    # to the subset of atoms to be included in B_damage calculations. The
    # subset of atoms retained is specified by the 'HETATM',
    # 'proteinOrNucleicAcid', 'addAtoms'  and 'removeAtoms' argument values
    # specified in INPUT.txt:

    import copy
    duplicate = copy.copy

    bdamatomList_1 = atomList
    bdamatomList_2 = duplicate(bdamatomList_1)

    # Removes hetatm if HETATM set to 'Remove' in INPUT.txt.
    for index, atm in enumerate(bdamatomList_1):
        if atm.lineID == 'HETATM':
            if HETATM is False:
                bdamatomList_2[index] = None

    for index, atm in enumerate(bdamatomList_1):
        if atm.lineID == 'ATOM':
            # Removes nucleic acid atoms if proteinOrNucleicAcid set to
            # 'Protein' in INPUT.txt.
            if protOrNA == 'PROTEIN':
                if len(str(atm.resiType)) != 3:
                    bdamatomList_2[index] = None

            # Removes protein atoms if proteinOrNucleicAcid set to
            # 'Nucleic Acid' in INPUT.txt.
            elif protOrNA in ['NUCLEICACID', 'NA']:
                if len(str(atm.resiType)) == 3:
                    bdamatomList_2[index] = None

    # Removes atoms whose number is in removeAtoms list.
    for index, atm in enumerate(bdamatomList_1):
        if str(atm.atomNum) in removeAtoms:
            bdamatomList_2[index] = None

    # Removes atoms whose atom type is in removeAtoms list.
    for index, atm in enumerate(bdamatomList_1):
        if str(atm.atomType) in removeAtoms:
            bdamatomList_2[index] = None

    # Removes atoms whose residue type is in removeAtoms list.
    for index, atm in enumerate(bdamatomList_1):
        if str(atm.resiType) in removeAtoms:
            bdamatomList_2[index] = None

    bdamatomList_2 = filter(None, bdamatomList_2)
    bdal2Append = bdamatomList_2.append

    fileOpen = open(fileName, 'r')
    for line in fileOpen.readlines():
        if str((line[0:6]).strip()) in ['ATOM', 'HETATM']:
            # Adds atoms whose number is in addAtoms list.
            if str((line[6:11]).strip()) in addAtoms:
                y = atom()
                y.lineID = str(line[0:6].strip())
                y.atomNum = int(line[6:11].strip())
                y.atomType = str(line[12:16].strip())
                y.resiType = str(line[17:20].strip())
                y.chainID = str(line[21:22].strip())
                y.resiNum = int(line[22:26].strip())
                y.xyzCoords = [[float(line[30:38].strip())],
                               [float(line[38:46].strip())],
                               [float(line[46:54].strip())]]
                y.occupancy = float(line[54:60].strip())
                y.bFactor = float(line[60:66].strip())
                y.atomID = str(line[76:78].strip())
                y.charge = str(line[78:80].strip())
                bdal2Append(y)

            # Adds atoms whose atom type is in addAtoms list.
            elif str(((line[12:16]).strip()).upper()) in addAtoms:
                y = atom()
                y.lineID = str(line[0:6].strip())
                y.atomNum = int(line[6:11].strip())
                y.atomType = str(line[12:16].strip())
                y.resiType = str(line[17:20].strip())
                y.chainID = str(line[21:22].strip())
                y.resiNum = int(line[22:26].strip())
                y.xyzCoords = [[float(line[30:38].strip())],
                               [float(line[38:46].strip())],
                               [float(line[46:54].strip())]]
                y.occupancy = float(line[54:60].strip())
                y.bFactor = float(line[60:66].strip())
                y.atomID = str(line[76:78].strip())
                y.charge = str(line[78:80].strip())
                bdal2Append(y)

            # Adds atoms whose residue type is in addAtoms list.
            elif str(((line[17:20]).strip()).upper()) in addAtoms:
                y = atom()
                y.lineID = str(line[0:6].strip())
                y.atomNum = int(line[6:11].strip())
                y.atomType = str(line[12:16].strip())
                y.resiType = str(line[17:20].strip())
                y.chainID = str(line[21:22].strip())
                y.resiNum = int(line[22:26].strip())
                y.xyzCoords = [[float(line[30:38].strip())],
                               [float(line[38:46].strip())],
                               [float(line[46:54].strip())]]
                y.occupancy = float(line[54:60].strip())
                y.bFactor = float(line[60:66].strip())
                y.atomID = str(line[76:78].strip())
                y.charge = str(line[78:80].strip())
                bdal2Append(y)

    fileOpen.close()

    # Filters list of atoms for B_damage analysis to remove multiple copies
    # of the same atom.
    bdamatomSet = set()
    bdamatomList = []
    for atm in bdamatomList_2:
        if atm.atomNum not in bdamatomSet:
            bdamatomList.append(atm)
            bdamatomSet.add(atm.atomNum)

    print 'Finished reading in atoms --> %d atoms found' % int(len(atomList))
    print 'in %s' % fileName

    bdamatomList = sorted(bdamatomList, key=lambda x: x.atomNum)
    return bdamatomList


def getUnitCellParams(fileName):
    # Parses in the unit cell parameters (the vectors a, b, c, and the angles
    # alpha, beta and gamma) from the input PDB file.

    import math

    fileOpen = open(fileName, 'r')
    for line in fileOpen.readlines():
        if('CRYST1' in str(line[0:6])):  # Unit cell parameters are stored in this line
            params = line.split()
            a = float(params[1])
            b = float(params[2])
            c = float(params[3])
            alpha = math.radians(float(params[4]))
            beta = math.radians(float(params[5]))
            gamma = math.radians(float(params[6]))
            break

    print 'Unit cell parameters successfully extracted'
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


def trimAtoms(atomList, params):
    # Removes all atoms with coordinates which lie outside of the trimmed
    # atoms box from the list of atoms in the 3x3 unit cell assembly.

    from atomCheck import isInXYZparams

    print 'Excluding the atoms that lie outside of the box'

    trimAtomList = []
    trimAppend = trimAtomList.append

    for atom in atomList:
        atomXYZ = atom.xyzCoords
        if isInXYZparams(atomXYZ, params):
            trimAppend(atom)

    print '%.0f atoms have been retained\n' % int(len(trimAtomList))
    return trimAtomList
