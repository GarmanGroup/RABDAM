
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


class atom(object):
    # A class to assign properties to atoms (specifically all properties
    # provided in the input PDB / mmCif file, plus values relating to the
    # BDamage calculation).

    def __init__(self, lineidentifier='', atomnum=0, atomtype='', conformer='',
                 resitype='', chainID='', residuenum=0, insertioncode = '',
                 xyz_coords=[], atomidentifier='', bfactor=0, occupancy=1,
                 charge='', packingdensity=0, avrg_bfactor=0, bdamage=0):
        self.lineID = lineidentifier
        self.atomNum = atomnum
        self.atomType = atomtype
        self.conformer = conformer
        self.resiType = resitype
        self.chainID = chainID
        self.resiNum = residuenum
        self.insCode = insertioncode
        self.xyzCoords = xyz_coords
        self.occupancy = occupancy
        self.bFactor = bfactor
        self.atomID = atomidentifier
        self.charge = charge
        self.pd = packingdensity
        self.avrg_bf = avrg_bfactor
        self.bd = bdamage


def download_pdb_and_mmcif(PDBcode, PDBdirectory, pathToPDB, pathToCif):
    # Downloads and saves PDB and mmCif files from the RSCB PDB website.

    import os
    import requests

    pdb_url = 'http://www.rcsb.org/pdb/files/%s.pdb' % PDBcode
    cif_url = 'http://www.rcsb.org/pdb/files/%s.cif' % PDBcode

    os.makedirs(PDBdirectory)
    print('\nDirectory %s created' % PDBdirectory)

    origPDB = requests.get(pdb_url)
    print('Downloaded PDB file from %s' % pdb_url)
    pdb_file = open(pathToPDB, 'w')
    pdb_file.write(origPDB.text)
    print('PDB file saved to %s' % pathToPDB)
    pdb_file.close()

    origCif = requests.get(cif_url)
    print('Downloaded mmCif file from %s' % cif_url)
    cif_file = open(pathToCif, 'w')
    cif_file.write(origCif.text)
    print('mmCif file saved to %s' % pathToCif)
    cif_file.close()


def copy_input(pathToInput, disk, newPathToInput, PDBdirectory):
    # Copies specified file to Logfiles directory.

    import os

    owd = os.getcwd()
    os.chdir('/')  # Changes directory to root directory
    os.chdir(disk)
    origInput = open(pathToInput, 'r')
    os.chdir(owd)
    os.makedirs(PDBdirectory)
    localFile = open(newPathToInput, 'w')
    localFile.write(origInput.read())
    print('%s copied to %s' % (pathToInput, newPathToInput))
    localFile.close()
    origInput.close()


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
            new_atom.insCode = line[26:27].strip()
            new_atom.xyzCoords = [[float(line[30:38].strip())],
                                  [float(line[38:46].strip())],
                                  [float(line[46:54].strip())]]
            new_atom.occupancy = float(line[54:60].strip())
            new_atom.bFactor = float(line[60:66].strip())
            new_atom.atomID = line[76:78].strip()
            new_atom.charge = line[78:80].strip()
            ucAtomList.append(new_atom)
    fileOpen.close()
    print('Finished reading in atoms --> %d atoms found in\n%s' % (
          int(len(ucAtomList)), fileName))
    return ucAtomList


def b_damage_atom_list(clean_au_list, HETATM, protOrNA, addAtoms,
                       removeAtoms):
    # Filters a copy of the clean_au_list to retain only the subset of atoms
    # to be included in the BDamage calculation, as specified by the 'HETATM',
    # 'proteinOrNucleicAcid', 'addAtoms'  and 'removeAtoms' argument values set
    # in the input file.

    import copy
    duplicate = copy.copy

    bdam_list_unfiltered = duplicate(clean_au_list) + [None]*len(clean_au_list)

    for index, atm in enumerate(clean_au_list):
        # Removes hetatm if HETATM set to 'Remove' in input file.
        if atm.lineID == 'HETATM':
            if HETATM is False:
                if atm.resiType not in ['MSE']:
                    bdam_list_unfiltered[index] = None

        elif atm.lineID == 'ATOM':
            # Removes nucleic acid atoms if proteinOrNucleicAcid set to
            # 'Protein' in input file.
            if protOrNA == 'protein':
                if len(atm.resiType) != 3:
                    bdam_list_unfiltered[index] = None
            # Removes protein atoms if proteinOrNucleicAcid set to
            # 'Nucleic Acid' / 'NA' in input file.
            elif protOrNA in ['nucleicacid', 'na']:
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

    print('Finished reading in atoms --> %d atoms found' % int(len(clean_au_list)))
    return bdamatomList
