
# RABDAM
# Copyright (C) 2020 Garman Group, University of Oxford

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
                 packingdensity=None, avrg_bfactor=None, bdamage=None):
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
               self.bd == other.bd


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
    print('mmCIF file saved to %s' % pathToCIF)
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
    print('%s copied to %s' % (pathToInput, newPathToInput))
    localFile.close()
    origInput.close()


def full_atom_list(fileName):
    """
    Parses in input PDB file (= PDBCUR output) and returns a list of all
    atoms in the file with their associated attributes as assigned in the
    'atom' class.
    """

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


def b_damage_atom_list(clean_au_list, seqres, HETATM, protOrNA, addAtoms,
                       removeAtoms):
    """
    Filters a copy of the clean_au_list to retain only the subset of atoms
    to be included in the BDamage calculation, as specified by the 'HETATM',
    'proteinOrNucleicAcid', 'addAtoms'  and 'removeAtoms' argument values set
    in the input file.
    """

    import copy

    if __name__ == 'Subroutines.parsePDB':
        from Subroutines.check_chem_components import nuc_acid_codes, amino_acid_codes
    else:
        from rabdam.Subroutines.check_chem_components import nuc_acid_codes, amino_acid_codes

    all_na_codes = nuc_acid_codes()
    all_aa_codes = amino_acid_codes()
    std_na_codes = ['A', 'C', 'G', 'I', 'U', 'DA', 'DC', 'DG', 'DI', 'DT',
                    'DU', 'N']
    std_aa_codes = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
                    'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
                    'THR', 'TRP', 'TYR', 'VAL', 'PYL', 'SEC', 'ASX', 'GLX',
                    'UNK']
    std_codes = std_na_codes + std_aa_codes

    bdam_list_unfiltered = copy.deepcopy(clean_au_list) + [None]*len(clean_au_list)

    for index, atm in enumerate(clean_au_list):
        # Removes hetatm not in macromolecule if HETATM set to 'Remove' in
        # input file.
        if atm.lineID == 'HETATM':
            if HETATM is False:
                if atm.resiType not in seqres:
                    bdam_list_unfiltered[index] = None
                elif (    atm.resiType in seqres
                      and any(x == atm.resiType for x in std_codes)
                ):
                    bdam_list_unfiltered[index] = None

        # Removes nucleic acid atoms if proteinOrNucleicAcid set to
        # 'Protein' in input file.
        if protOrNA == 'protein':
            if atm.resiType in seqres and not atm.resiType in all_aa_codes:
                bdam_list_unfiltered[index] = None
        # Removes protein atoms if proteinOrNucleicAcid set to
        # 'Nucleic Acid' / 'NA' in input file.
        elif protOrNA in ['nucleicacid', 'na']:
            if atm.resiType in seqres and not atm.resiType in all_na_codes:
                bdam_list_unfiltered[index] = None

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
