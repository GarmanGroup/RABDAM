
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


def convert_array_to_atom_list(array, atom_id_list, atom_list):
    # Converts numpy array of xyz coordinates into a list of Atom() objects,
    # for use in writing a PDB file

    import copy
    duplicate = copy.copy
    import numpy as np

    pdb_atom_list = []
    for index, atom_id in enumerate(atom_id_list):
        for atom in atom_list:
            if atom.atomNum == atom_id:
                new_atom = duplicate(atom)
                newXYZ = np.array([[array[index][0]],
                                   [array[index][1]],
                                   [array[index][2]]])
                new_atom.xyzCoords = newXYZ
                pdb_atom_list.append(new_atom)
                break

    return pdb_atom_list


def makePDB(header_lines, atomList, footer_lines, newPDBfilename, Bfac):
    # Writes a PDB file containing a complete set of atom information for all
    # atoms in 'atomList', plus header and footer information.

    import numpy as np

    newPDBfile = open(newPDBfilename, 'a')

    for line in header_lines:
        newPDBfile.write(line)

    for index, atm in enumerate(atomList):
        a = atm.lineID
        b = atm.atomNum
        c = atm.atomType
        d = atm.conformer
        e = atm.resiType
        f = atm.chainID
        g = atm.resiNum
        h = atm.insCode
        i = atm.xyzCoords[0][0]
        j = atm.xyzCoords[1][0]
        k = atm.xyzCoords[2][0]
        l = atm.occupancy
        if Bfac == 'Bfactor':
            m = atm.bFactor
        elif Bfac == 'BDamage':
            m = np.log(atm.bd)
        n = atm.atomID
        o = atm.charge
        # Atom properties are appropriately ordered and spaced, and reported
        # to the expected number of significant figures, for the PDB file
        # format. Note that atomType for some metal ions will not follow
        # standard PDB file format, but this will not affect the running of
        # RABDAM (nor most other programs that the user might want to load the
        # PDB file into, such as PyMol, Chimera, CCP4MG, WinCoot, etc.)
        newLine = '%-6s%5d  %-3s%1s%3s%2s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n' % (
            a, b, c, d, e, f, g, h, i, j, k, l, m, n, o
            )
        newPDBfile.write(newLine)

        # Inserts TER cards
        if index != (len(atomList) - 1):
            if (
                atm.chainID != atomList[index+1].chainID
                and (atm.resiNum+1) != atomList[index+1].resiNum
                and not (atm.lineID == 'HETATM' and atomList[index+1].lineID == 'HETATM')
                ):
                newPDBfile.write('TER\n')
        else:
            if atm.lineID != 'HETATM':
                newPDBfile.write('TER\n')

    for line in footer_lines:
        newPDBfile.write(line)

    print('New PDB file saved to %s' % newPDBfilename)
    newPDBfile.close()


def writeDataFrame(bdamAtomList):
    # Returns a DataFrame containing a complete set of atom information
    # (including both that provided in the input PDB file and also the
    # BDamage values calculated by RABDAM) for all atoms considered for
    # BDamage analysis.

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
        ATMNAME[index] = atm.atomType
        CONFORMER[index] = atm.conformer
        RESNAME[index] = atm.resiType
        CHAIN[index] = atm.chainID
        RESNUM[index] = atm.resiNum
        INSCODE[index] = atm.insCode
        XPOS[index] = atm.xyzCoords[0][0]
        YPOS[index] = atm.xyzCoords[1][0]
        ZPOS[index] = atm.xyzCoords[2][0]
        OCC[index] = atm.occupancy
        BFAC[index] = atm.bFactor
        ELEMENT[index] = atm.atomID
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

    df_list_dict_copy = copy.copy(df_list_dict)
    for key in df_list_dict_copy:
        if set(df_list_dict_copy[key]) == {''}:
            df_list_dict[key] = ['?']*len(bdamAtomList)
        elif (len(set(df_list_dict_copy[key])) > 1
            and '' in set(df_list_dict_copy[key])
            ):
            cif_list = []
            for atm in df_list_dict_copy[key]:
                if atm == '':
                    cif_list.append('.')
                else:
                    cif_list.append(atm)
            df_list_dict[key] = cif_list

    # Generates dictionary of column widths for output cif file
    cif_column_widths = {'REC': len(max(df_list_dict['REC'], key=len)),
                         'ATMNUM': len(max([str(value) for value in df_list_dict['ATMNUM']], key=len)),
                         'ATMNAME': len(max(df_list_dict['ATMNAME'], key=len)),
                         'CONFORMER': len(max(df_list_dict['CONFORMER'], key=len)),
                         'RESNAME': len(max(df_list_dict['RESNAME'], key=len)),
                         'CHAIN': len(max(df_list_dict['CHAIN'], key=len)),
                         'RESNUM': len(max([str(value) for value in df_list_dict['RESNUM']], key=len)),
                         'INSCODE': len(max(df_list_dict['INSCODE'], key=len)),
                         'XPOS': len(max([str(value) for value in df_list_dict['XPOS']], key=len)),
                         'YPOS': len(max([str(value) for value in df_list_dict['YPOS']], key=len)),
                         'ZPOS': len(max([str(value) for value in df_list_dict['ZPOS']], key=len)),
                         'OCC': len(max([str(value) for value in df_list_dict['OCC']], key=len)),
                         'BFAC': len(max([str(value) for value in df_list_dict['BFAC']], key=len)),
                         'BDAM': len(max(['{0:.2f}'.format(value) for value in df_list_dict['BDAM']], key=len)),
                         'ELEMENT': len(max(df_list_dict['ELEMENT'], key=len)),
                         'CHARGE': len(max(df_list_dict['CHARGE'], key=len))}

    # Lists are concatenated into the colummns of a DataFrame.
    df = pd.DataFrame(df_list_dict)

    # DataFrame columns are ordered.
    cols = df.columns.tolist()
    cols = ([cols[12]] + [cols[1]] + [cols[0]] + [cols[7]] + [cols[13]]
            + [cols[5]] + [cols[14]] + [cols[9]] + [cols[15]] + [cols[16]]
            + [cols[17]] + [cols[10]] + [cols[4]] + [cols[8]] + [cols[6]]
            + [cols[11]] + [cols[2]] + [cols[3]])
    df = df[cols]

    return df, cif_column_widths
