
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


def convert_array_to_atom_list(array, all_atoms_list, sub_atoms_list):
    """
    Converts numpy array of xyz coordinates into a list of Atom() objects,
    for use in writing a PDB file
    """

    import copy
    import numpy as np

    pdb_atom_list = []
    for index, atom_id in enumerate(all_atoms_list):
        for atom in sub_atoms_list:
            if atom.atomNum == atom_id:
                new_atom = copy.deepcopy(atom)
                newXYZ = np.array([[array[index][0]],
                                   [array[index][1]],
                                   [array[index][2]]])
                new_atom.xyzCoords = newXYZ
                pdb_atom_list.append(new_atom)
                break

    return pdb_atom_list


def makePDB(header_lines, atomList, footer_lines, seqres, newPDBfilename, Bfac):
    """
    Writes a PDB file containing a complete set of atom information for all
    atoms in 'atomList', plus header and footer information.
    """

    import numpy as np

    exit = False
    newPDBfile = open(newPDBfilename, 'w')

    for line in header_lines:
        if not line.endswith('\n'):
            line += '\n'
        newPDBfile.write(line)

    for index, atm in enumerate(atomList):
        a = atm.lineID.ljust(6)
        if len(a) != 6:
            exit = True
        b = str(atm.atomNum).rjust(5)[-5:]
        c = atm.atomType.ljust(3)
        if len(c) != 3:
            exit = True
        d = atm.conformer.ljust(1)
        if len(d) != 1:
            exit = True
        e = atm.resiType.ljust(3)
        if len(e) != 3:
            exit = True
        f = atm.chainID.ljust(1)
        if len(f) != 1:
            exit = True
        g = str(atm.resiNum).rjust(4)[-4:]
        h = atm.insCode.ljust(1)
        if len(h) != 1:
            exit  = True
        i = '%8.3f' % atm.xyzCoords[0][0]
        if len(i) != 8:
            exit = True
        j = '%8.3f' % atm.xyzCoords[1][0]
        if len(j) != 8:
            exit = True
        k = '%8.3f' % atm.xyzCoords[2][0]
        if len(k) != 8:
            exit = True
        l = '%6.2f' % atm.occupancy
        if len(l) != 6:
            exit = True
        if Bfac.lower() == 'bfactor':
            m = '%6.2f' % atm.bFactor
            if len(m) != 6:
                exit = True
        elif Bfac.lower() == 'bdamage':
            m = '%6.2f' % np.log(atm.bd)
            if len(m) != 6:
                exit = True
        n = atm.element.rjust(2)
        if len(n) != 2:
            exit = True
        o = atm.charge.rjust(2)
        if len(o) != 2:
            exit = True

        if exit is False:
            # Atom properties are appropriately ordered and spaced, and reported
            # to the expected number of significant figures, for the PDB file
            # format. Note that atomType for some metal ions will not follow
            # standard PDB file format, but this will not affect the running of
            # RABDAM (nor most other programs that the user might want to load
            # the PDB file into, such as PyMol, Chimera, CCP4MG, WinCoot, etc.)
            newLine = ''.join([a, b, '  ', c, d, e, ' ', f, g, h, '   ', i, j,
                               k, l, m, '          ', n, o, '\n'])
            newPDBfile.write(newLine)

            # Inserts TER cards
            if index != (len(atomList) - 1):
                if (
                        (atm.chainID != atomList[index+1].chainID)
                    and (atm.resiType in seqres)
                ):
                    newPDBfile.write('TER'.ljust(80) + '\n')
            else:
                if atm.resiType in seqres:
                    newPDBfile.write('TER'.ljust(80) + '\n')

    for line in footer_lines:
        if not line.endswith('\n'):
            line += '\n'
        newPDBfile.write(line)

    print('New PDB file saved to %s' % newPDBfilename)
    newPDBfile.close()

    return exit


def writeDataFrame(bdamAtomList):
    """
    Returns a DataFrame containing a complete set of atom information
    (including both that provided in the input PDB file and also the BDamage
    values calculated by RABDAM) for all atoms considered for BDamage analysis.
    """

    import copy
    import pandas as pd

    # Initialises a list for each atom property considered.
    REC = [None]*len(bdamAtomList)
    ATMNUM = [None]*len(bdamAtomList)
    ATMNAME = [None]*len(bdamAtomList)
    CONFORMER = [None]*len(bdamAtomList)
    RESNAME = [None]*len(bdamAtomList)
    CHAIN = [None]*len(bdamAtomList)
    RESNUM = [None]*len(bdamAtomList)
    INSCODE = [None]*len(bdamAtomList)
    XPOS = [None]*len(bdamAtomList)
    YPOS = [None]*len(bdamAtomList)
    ZPOS = [None]*len(bdamAtomList)
    OCC = [None]*len(bdamAtomList)
    BFAC = [None]*len(bdamAtomList)
    ELEMENT = [None]*len(bdamAtomList)
    CHARGE = [None]*len(bdamAtomList)
    PD = [None]*len(bdamAtomList)
    AVRG_BF = [None]*len(bdamAtomList)
    BDAM = [None]*len(bdamAtomList)

    # Lists are filled with the relevant values of the properties associated
    # with each of the atoms considered for BDamage analysis.
    for index, atm in enumerate(bdamAtomList):
        REC[index] = atm.lineID
        ATMNUM[index] = atm.atomNum
        ATMNAME[index] = atm.origAtomType
        CONFORMER[index] = atm.conformer
        RESNAME[index] = atm.origResiType
        CHAIN[index] = atm.origChainID
        RESNUM[index] = atm.origResiNum
        INSCODE[index] = atm.insCode
        XPOS[index] = atm.xyzCoords[0][0]
        YPOS[index] = atm.xyzCoords[1][0]
        ZPOS[index] = atm.xyzCoords[2][0]
        OCC[index] = atm.occupancy
        BFAC[index] = atm.bFactor
        ELEMENT[index] = atm.element
        CHARGE[index] = atm.charge
        PD[index] = atm.pd
        AVRG_BF[index] = atm.avrg_bf
        BDAM[index] = atm.bd

    # Generates dictionary of DataFrame columns
    df_list_dict = {'REC': REC,
                    'ATMNUM': ATMNUM,
                    'ATMNAME': ATMNAME,
                    'CONFORMER': CONFORMER,
                    'RESNAME': RESNAME,
                    'CHAIN': CHAIN,
                    'RESNUM': RESNUM,
                    'INSCODE': INSCODE,
                    'XPOS': XPOS,
                    'YPOS': YPOS,
                    'ZPOS': ZPOS,
                    'OCC': OCC,
                    'BFAC': BFAC,
                    'ELEMENT': ELEMENT,
                    'CHARGE': CHARGE,
                    'PD': PD,
                    'AVRG_BF': AVRG_BF,
                    'BDAM': BDAM}

    # Ensures output dataframe is in cif format
    df_list_dict_copy = copy.copy(df_list_dict)
    for key in df_list_dict_copy:
        if set(df_list_dict_copy[key]) == {''}:
            df_list_dict[key] = ['?']*len(bdamAtomList)
        elif (
                len(set(df_list_dict_copy[key])) > 1
            and '' in set(df_list_dict_copy[key])
        ):
            cif_list = []
            for atm in df_list_dict_copy[key]:
                if atm == '':
                    cif_list.append('.')
                else:
                    cif_list.append(atm)
            df_list_dict[key] = cif_list

    # Lists are concatenated into the colummns of a DataFrame.
    df = pd.DataFrame(df_list_dict)

    # DataFrame columns are ordered.
    df = df[['REC', 'ATMNUM', 'ATMNAME', 'CONFORMER', 'RESNAME', 'CHAIN',
             'RESNUM', 'INSCODE', 'XPOS', 'YPOS', 'ZPOS', 'OCC', 'BFAC',
             'ELEMENT', 'CHARGE', 'PD', 'AVRG_BF', 'BDAM']]

    return df
