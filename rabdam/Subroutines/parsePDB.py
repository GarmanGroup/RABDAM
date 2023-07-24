
# RABDAM
# Copyright (C) 2023 Garman Group, University of Oxford

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
    """
    A class to assign properties to atoms (specifically all properties
    provided in the input PDB / mmCif file, plus values relating to the
    BDamage calculation).
    """

    def __init__(self, lineidentifier=None, atomnum=None, atomtype=None,
                 conformer=None, resitype=None, chainID=None, resinum=None,
                 insertioncode=None, xyz_coords=None, element=None,
                 bfactor=None, occupancy=None, charge=None, orig_resinum=None,
                 orig_resitype=None, orig_chainID=None, orig_atomtype=None,
                 packingdensity=None, avrg_bfactor=None, bdamage=None,
                 protein=None, na=None, chain_len=None):
        self.lineID = lineidentifier
        self.atomNum = atomnum
        self.atomType = atomtype
        self.conformer = conformer
        self.resiType = resitype
        self.chainID = chainID
        self.resiNum = resinum
        self.insCode = insertioncode
        self.xyzCoords = xyz_coords
        self.occupancy = occupancy
        self.bFactor = bfactor
        self.element = element
        self.charge = charge
        self.origResiNum = orig_resinum
        self.origResiType = orig_resitype
        self.origChainID = orig_chainID
        self.origAtomType = orig_atomtype
        self.pd = packingdensity
        self.avrg_bf = avrg_bfactor
        self.bd = bdamage
        self.protein = protein
        self.na = na
        self.chain_len = chain_len

    def __eq__(self, other):
        """
        Allows atom objects to be compared during unit testing
        """

        return self.lineID == other.lineID and \
               self.atomNum == other.atomNum and \
               self.atomType == other.atomType and \
               self.conformer == other.conformer and \
               self.resiType == other.resiType and \
               self.chainID == other.chainID and \
               self.resiNum == other.resiNum and \
               self.insCode == other.insCode and \
               self.xyzCoords == other.xyzCoords and \
               self.occupancy == other.occupancy and \
               self.bFactor == other.bFactor and \
               self.element == other.element and \
               self.charge == other.charge and \
               self.origResiNum == other.origResiNum and \
               self.origResiType == other.origResiType and \
               self.origChainID == other.origChainID and \
               self.origAtomType == other.origAtomType and \
               self.pd == other.pd and \
               self.avrg_bf == other.avrg_bf and \
               self.bd == other.bd and \
               self.protein == other.protein and\
               self.na == other.na and \
               self.chain_len == other.chain_len


def download_mmcif(PDBcode, PDBdirectory, pathToCIF):
    """
    Downloads and saves mmCif file from the RSCB PDB website.
    """

    import os
    import requests

    # Checks whether accession code exists - if not, exit program
    # with error message
    exit = False
    mmcif_url = 'https://files.rcsb.org/view/%s.cif' % PDBcode
    header = requests.get(mmcif_url)
    if header.status_code != 200:
        print('\n\nERROR: Failed to download %s mmCIF file with accession '
              'code %s:\nCheck that a structure with this accession code '
              'exists.' % (mmcif_url[-4:], PDBcode))
        exit = True

    # Downloads and saves mmCIF file
    os.makedirs(PDBdirectory)
    print('\nDirectory %s created' % PDBdirectory)

    origCIF = requests.get(mmcif_url)
    print('Downloaded mmCIF file from %s' % mmcif_url)
    cif_file = open(pathToCIF, 'w')
    cif_file.write(origCIF.text)
    print('mmCIF file saved to %s/%s' % (PDBdirectory, pathToCIF))
    cif_file.close()

    return exit


def copy_input(pathToInput, disk, newPathToInput, PDBdirectory):
    """
    Copies specified file to Logfiles directory.
    """

    import os

    owd = os.getcwd()
    os.chdir('/')  # Changes directory to root directory
    os.chdir(disk)
    origInput = open(pathToInput, 'r')
    os.chdir(owd)
    os.makedirs(PDBdirectory)
    localFile = open(newPathToInput, 'w')
    localFile.write(origInput.read())
    print('%s copied to %s/%s' % (pathToInput, owd, newPathToInput))
    localFile.close()
    origInput.close()


def b_damage_atom_list(clean_au_list, HETATM, protOrNA, addAtoms, removeAtoms):
    """
    Filters a copy of the clean_au_list to retain only the subset of atoms
    to be included in the BDamage calculation, as specified by the 'HETATM',
    'proteinOrNucleicAcid', 'addAtoms'  and 'removeAtoms' argument values set
    in the input file.
    """

    import copy

    bdam_list_unfiltered = copy.deepcopy(clean_au_list) + [None]*len(clean_au_list)

    for index, atm in enumerate(clean_au_list):
        # Removes hetatm not in macromolecule if HETATM set to 'Remove' in
        # input file.
        if atm.protein is False and atm.na is False:
            if HETATM is False:
                bdam_list_unfiltered[index] = None

        # Removes nucleic acid atoms if proteinOrNucleicAcid set to
        # 'Protein' in input file.
        if protOrNA == 'protein':
            if atm.protein is False and atm.na is True:
                bdam_list_unfiltered[index] = None
        # Removes protein atoms if proteinOrNucleicAcid set to
        # 'Nucleic Acid' / 'NA' in input file.
        elif protOrNA in ['nucleicacid', 'na']:
            if atm.na is False and atm.protein is True:
                bdam_list_unfiltered[index] = None
        # Otherwise keeps all protein and NA atoms
        elif protOrNA == 'proteinna':
            pass

        # Removes atoms whose number is in removeAtoms list.
        if str(atm.atomNum) in removeAtoms:
            bdam_list_unfiltered[index] = None
        # Removes atoms whose residue type is in removeAtoms list.
        elif atm.resiType in removeAtoms:
            bdam_list_unfiltered[index] = None

        # Adds atoms whose number is in addAtoms list.
        if str(atm.atomNum) in addAtoms:
            bdam_list_unfiltered[len(clean_au_list)+index] = copy.deepcopy(atm)
        # Adds atoms whose residue type is in addAtoms list.
        elif atm.resiType in addAtoms:
            bdam_list_unfiltered[len(clean_au_list)+index] = copy.deepcopy(atm)

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
