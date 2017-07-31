

class atom(object):
    # A class to assign properties to atoms (specifically all properties
    # provided in the input PDB file, plus values relating to the BDamage
    # calculation).

    def __init__(self, lineidentifier='', atomnum=0, atomtype='', conformer='',
                 resitype='', chainID='', residuenum=0, xyz_coords=[],
                 atomidentifier='', bfactor=0, occupancy=1, charge='',
                 packingdensity=0, avrg_bfactor=0, bdamage=0):
        self.lineID = lineidentifier
        self.atomNum = atomnum
        self.atomType = atomtype
        self.conformer = conformer
        self.resiType = resitype
        self.chainID = chainID
        self.resiNum = residuenum
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
    import requests

    urlText = 'http://www.rcsb.org/pdb/files/%s.pdb' % PDBcode
    os.makedirs(PDBdirectory)
    print 'Directory %s created' % PDBdirectory
    origPDB = requests.get(urlText)
    print 'Downloaded PDB file from %s' % urlText
    localFile = open(pathToPDB, 'w')
    localFile.write(origPDB.text)
    print 'PDB file saved to %s' % pathToPDB
    localFile.close()


def copyPDB(pathToPDB, disk, newPathToPDB, PDBdirectory):
    # Copies .pdb / .txt file from any provided file path to Logfiles
    # directory.

    import os

    owd = os.getcwd()
    os.chdir('/')  # Changes directory to root directory
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
    # 'atom' class.

    ucAtomList = []
    fileOpen = open(fileName, 'r')
    for line in fileOpen.readlines():
        if line[0:6].strip() in ['ATOM', 'HETATM']:
            new_atom = atom()
            new_atom.lineID = line[0:6].strip()
            new_atom.atomNum = int(line[6:11].strip())
            new_atom.atomType = line[12:16].strip()
            new_atom.conformer = line[16:17].strip()
            new_atom.resiType = line[17:20].strip()
            new_atom.chainID = line[21:22].strip()
            new_atom.resiNum = int(line[22:26].strip())
            new_atom.xyzCoords = [[float(line[30:38].strip())],
                                  [float(line[38:46].strip())],
                                  [float(line[46:54].strip())]]
            new_atom.occupancy = float(line[54:60].strip())
            new_atom.bFactor = float(line[60:66].strip())
            new_atom.atomID = line[76:78].strip()
            new_atom.charge = line[78:80].strip()
            ucAtomList.append(new_atom)
    fileOpen.close()
    print 'Finished reading in atoms --> %d atoms found' % int(len(ucAtomList))
    print 'in %s' % fileName
    return ucAtomList


def b_damage_atom_list(clean_au_list, HETATM, protOrNA, addAtoms,
                       removeAtoms):
    # Filters a copy of the clean_au_list to retain only the subset of atoms
    # to be included in the BDamage calculation, as specified by the 'HETATM',
    # 'proteinOrNucleicAcid', 'addAtoms'  and 'removeAtoms' argument values set
    # in the input file.

    import copy
    duplicate = copy.copy

    bdam_list_unfiltered = duplicate(clean_au_list)
    bdam_list_unfiltered = bdam_list_unfiltered + [None]*len(bdam_list_unfiltered)

    for index, atm in enumerate(clean_au_list):
        # Removes hetatm if HETATM set to 'Remove' in input file.
        if atm.lineID == 'HETATM':
            if HETATM is False:
                bdam_list_unfiltered[index] = None

        elif atm.lineID == 'ATOM':
            # Removes nucleic acid atoms if proteinOrNucleicAcid set to
            # 'Protein' in input file.
            if protOrNA == 'PROTEIN':
                if len(atm.resiType) != 3:
                    bdam_list_unfiltered[index] = None
            # Removes protein atoms if proteinOrNucleicAcid set to
            # 'Nucleic Acid' / 'NA' in input file.
            elif protOrNA in ['NUCLEICACID', 'NA']:
                if len(atm.resiType) == 3:
                    bdam_list_unfiltered[index] = None

        # Removes atoms whose number is in removeAtoms list.
        if str(atm.atomNum) in removeAtoms:
            bdam_list_unfiltered[index] = None
        # Removes atoms whose residue type is in removeAtoms list.
        elif atm.resiType in removeAtoms:
            bdam_list_unfiltered[index] = None

        # Adds atoms whose number is in addAtoms list.
        if str(atm.atomNum) in addAtoms:
            bdam_list_unfiltered[len(clean_au_list)+index] = atm
        # Adds atoms whose residue type is in addAtoms list.
        elif atm.resiType in addAtoms:
            bdam_list_unfiltered[len(clean_au_list)+index] = atm

    bdam_list_filtered = filter(None, bdam_list_unfiltered)
    bdam_list_filtered = sorted(bdam_list_filtered, key=lambda x: x.atomNum)

    # Filters list of atoms for BDamage analysis to remove multiple copies
    # of the same atom.
    bdamatomSet = set()
    bdamatomList = []
    for atm in bdam_list_filtered:
        if atm.atomNum not in bdamatomSet:
            bdamatomList.append(atm)
            bdamatomSet.add(atm.atomNum)

    print 'Finished reading in atoms --> %d atoms found' % int(len(clean_au_list))
    return bdamatomList