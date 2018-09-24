
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


def convert_cif_to_pdb(pathToInput, convert_cif):
    # If RABDAM is provided with an input mmCIF file, converts it into PDB file
    # format.
    import os
    import sys

    if __name__ == 'Subroutines.PDBCUR':
        from Subroutines.parsePDB import atom
        from Subroutines.makeDataFrame import makePDB
    else:
        from rabdam.Subroutines.parsePDB import atom
        from rabdam.Subroutines.makeDataFrame import makePDB

    exit = False

    input_cif = open('%s' % pathToInput, 'r')
    input_cif_lines = ''
    for line in input_cif:
        line = line.replace('\n', '')
        line = line.replace('\r', '')
        line = line.strip()
        input_cif_lines = input_cif_lines + line + '\\'
    input_cif.close()

    input_cif_lines = input_cif_lines.split('#\\')
    ATM_rec = ''
    aniso_rec = ''
    disulfide_bonds = ''
    space_group = ''
    cif_header_lines = ''
    cif_footer_lines = '#\\'
    header = True
    footer = False
    for subsection in input_cif_lines:
        if '_atom_site.group_PDB' in subsection:
            ATM_rec = subsection
            header = False
        elif '_atom_site_anisotrop.id' in subsection:
            aniso_rec = subsection
            footer = True
        elif 'disulf1' in subsection:
            disulfide_bonds = subsection
        elif ('_cell.length_a' in subsection
             or '_symmetry.space_group_name_H-M' in subsection
             ):
            space_group = space_group + subsection

        if header is True:
            cif_header_lines = cif_header_lines + subsection + '#\\'
        elif footer is True and not '_atom_site_anisotrop.id' in subsection:
            cif_footer_lines = cif_footer_lines + subsection + '#\\'

    aniso_rec = aniso_rec.split('\\')
    aniso_rec.remove('')
    cif_header_lines = cif_header_lines.split('\\')
    cif_header_lines.remove('')
    cif_footer_lines = cif_footer_lines.split('\\')
    cif_footer_lines.remove('')

    # Constructs CRYST1 line of PDB file
    space_group = space_group.split('\\')
    a = ''
    b = ''
    c = ''
    alpha = ''
    beta = ''
    gamma = ''
    sGroup = ''
    z = ''
    for line in space_group:
        if line.startswith('_cell.length_a '):
            a = str(round(float(line.split()[-1].strip()), 3))
            if len(a) > 8:
                decimal_places = 7-len(a.split('.')[0])
                a = str(round(float(a), decimal_places))
        elif line.startswith('_cell.length_b '):
            b = str(round(float(line.split()[-1].strip()), 3))
            if len(b) > 8:
                decimal_places = len(b.split('.')[0])
                b = str(round(float(b), 7-decimal_places))
        elif line.startswith('_cell.length_c '):
            c = str(round(float(line.split()[-1].strip()), 3))
            if len(c) > 8:
                decimal_places = 7-len(c.split('.')[0])
                c = str(round(float(c), decimal_places))
        elif line.startswith('_cell.angle_alpha '):
            alpha = str(round(float(line.split()[-1].strip()), 2))
            if len(alpha) > 6:
                decimal_places = 5-len(alpha.split('.')[0])
                alpha = str(round(float(alpha), decimal_places))
        elif line.startswith('_cell.angle_beta '):
            beta = str(round(float(line.split()[-1].strip()), 2))
            if len(beta) > 6:
                decimal_places = 5-len(beta.split('.')[0])
                beta = str(round(float(beta), decimal_places))
        elif line.startswith('_cell.angle_gamma '):
            gamma = str(round(float(line.split()[-1].strip()), 2))
            if len(gamma) > 6:
                decimal_places = 5-len(gamma.split('.')[0])
                gamma = str(round(float(gamma), decimal_places))
        elif line.startswith('_symmetry.space_group_name_H-M '):
            sGroup = line.replace('_symmetry.space_group_name_H-M', '')
            sGroup = sGroup.replace("'", '')
            sGroup = sGroup.strip()
        elif line.startswith('_cell.Z_PDB '):
            z = line.split()[-1].strip()
    if any(x == '' for x in [a, b, c, alpha, beta, gamma, sGroup, z]):
        exit = True
        print('\n\nERROR: Failed to construct CRYST1 line')
    cryst1_line = ''.join(['CRYST1', a.rjust(9), b.rjust(9), c.rjust(9),
                           alpha.rjust(7), beta.rjust(7), gamma.rjust(7), ' ',
                           sGroup.ljust(11), z.rjust(4), '\n'])

    # Constructs SSBOND lines of PDB file
    disulfide_bonds = disulfide_bonds.split('\\')
    disulf_bonds = [line for line in disulfide_bonds if line.startswith('disulf')]
    disulf_column_labels = [line.replace('_struct_conn.', '') for line in
                            disulfide_bonds if line.startswith('_struct_conn.')]
    disulf_lines = ''
    for line in disulf_bonds:
        values = line.split()

        res_name_1 = values[disulf_column_labels.index('ptnr1_label_comp_id')]
        chain_id_1 = values[disulf_column_labels.index('ptnr1_label_asym_id')]
        res_num_1 = values[disulf_column_labels.index('ptnr1_label_seq_id')]
        ins_code_1 = values[disulf_column_labels.index('pdbx_ptnr1_PDB_ins_code')]
        res_name_2 = values[disulf_column_labels.index('ptnr2_label_comp_id')]
        chain_id_2 = values[disulf_column_labels.index('ptnr2_label_asym_id')]
        res_num_2 = values[disulf_column_labels.index('ptnr2_label_seq_id')]
        ins_code_2 = values[disulf_column_labels.index('pdbx_ptnr2_PDB_ins_code')]

        disulf_line = ''.join(['SSBOND', ''.rjust(5), res_name_1.rjust(3), ' ',
                               chain_id_1.rjust(1), ' ', res_num_1.rjust(4),
                               ins_code_1.rjust(1), '   ', res_name_2.rjust(3),
                               ' ', chain_id_2.rjust(1), ' ',
                               res_num_2.rjust(4), ins_code_2.rjust(1), '\n'])
        disulf_lines = disulf_lines + disulf_line

    # Constructs ATOM / HETATM record lines of PDB file
    if ATM_rec == '':
        exit = True
        print('\n\nERROR: Failed to extract ATOM / HETATM records from input '
              '.cif file\n- check that this file is consistent with the '
              'standard cif file format')
    else:
        if convert_cif is True:
            ATM_rec = ATM_rec.split('\\')
            cif_column_labels = [line.replace('_atom_site.', '') for line in
                                 ATM_rec if line.startswith('_atom_site.')]
            rows = [line for line in ATM_rec if line.startswith(('ATOM', 'HETATM'))]

            cif_atom_list = []
            multi_model = False
            for line in rows:
                values = line.split()

                if int(values[cif_column_labels.index('pdbx_PDB_model_num')]) > 1:
                    exit = True
                    print('\n\nERROR: More than one model present in input '
                          'mmCIF file.\nPlease submit an mmCIF file containing '
                          'a single model for BDamage analysis.\n')

                try:
                    input_atom = atom()
                    input_atom.lineID = values[cif_column_labels.index('group_PDB')]
                    input_atom.atomNum = int(values[cif_column_labels.index('id')])
                    input_atom.atomType = values[cif_column_labels.index('auth_atom_id')].strip('"')
                    input_atom.conformer = values[cif_column_labels.index('label_alt_id')]
                    input_atom.resiType = values[cif_column_labels.index('auth_comp_id')]
                    input_atom.chainID = values[cif_column_labels.index('auth_asym_id')]
                    input_atom.resiNum = int(values[cif_column_labels.index('auth_seq_id')])
                    input_atom.insCode = values[cif_column_labels.index('pdbx_PDB_ins_code')]
                    input_atom.xyzCoords = [[float(values[cif_column_labels.index('Cartn_x')])],
                                            [float(values[cif_column_labels.index('Cartn_y')])],
                                            [float(values[cif_column_labels.index('Cartn_z')])]]
                    input_atom.occupancy = float(values[cif_column_labels.index('occupancy')])
                    input_atom.bFactor = float(values[cif_column_labels.index('B_iso_or_equiv')])
                    input_atom.element = values[cif_column_labels.index('type_symbol')]
                    input_atom.charge = values[cif_column_labels.index('pdbx_formal_charge')]

                    cif_atom_list.append(input_atom)

                except ValueError:
                    exit = True
                    print('\n\nERROR: mmCIF file _atom_site labels do not '
                          'follow the expected format\n.')


            os.remove(pathToInput)
            pathToInput = pathToInput.replace('.cif', '.pdb')
            pdb_header_lines = cryst1_line + disulf_lines
            pdb_footer_lines = ''
            makePDB(pdb_header_lines, cif_atom_list, pdb_footer_lines,
                    pathToInput, 'Bfactor')

    return (pathToInput, aniso_rec, cif_header_lines, cif_footer_lines, exit)


def clean_pdb_file(pathToInput, PDBdirectory, pdb_file_path):
    # Filters the input PDB file ATOM / HETATM records to remove anisotropic
    # Bfactor records, hydrogen atoms and 0 occupancy atoms, as well as
    # retaining only the most probable alternate conformers (in the case
    # where multiple alternate conformers of a particular atom have the same
    # occupancy, (only) that which is first listed in the PDB file is
    # retained). The header and footer lines of the input PDB file are stored
    # in lists for use in writing output PDB files later in the program.

    import sys
    import shutil
    import math
    import numpy as np
    import pandas as pd

    if __name__ == 'Subroutines.PDBCUR':
        from Subroutines.parsePDB import atom
        from Subroutines.makeDataFrame import makePDB
    else:
        from rabdam.Subroutines.parsePDB import atom
        from rabdam.Subroutines.makeDataFrame import makePDB

    # Checks that only a single model is present in the input PDB file, and
    # that any disulphide bonds have been refined with 100% occupancy. Then
    # removes anisotropic Bfactors, hydrogen atoms and 0 occupancy atoms.
    exit = False
    pause = False

    filtered_pdb_lines = []
    header_lines = []
    footer_lines = []
    unit_cell_params = []

    orig_pdb = open('%s' % pathToInput, 'r')
    orig_pdb_lines = orig_pdb.readlines()
    header = True
    footer = False

    disulfide_bonds = []
    seqres = []
    atom_ids = []
    conformers = [[] for i in range(len(orig_pdb_lines))]
    occupancies = [[] for i in range(len(orig_pdb_lines))]

    for line in orig_pdb_lines:
        # Checks for single model
        if line.replace(' ', '').startswith('MODEL2'):
            exit = True
            print('\n\nERROR: More than one model present in input PDB file.\n'
                  'Please submit a PDB file containing a single model for '
                  'BDamage analysis.\n'
                  'Terminating RABDAM run.\n')

        # Extracts list of macromolecular residue ids
        elif line.startswith('SEQRES'):
            chain_id = line[11:12]
            res = ['{}_{}'.format(chain_id, resid) for resid in
                   line[19:].split() if not resid in ['', ' ', '\n']]
            seqres += res

        # Extracts ids of sulfur atoms involved in disulfide bonds
        elif line.startswith('SSBOND'):
            disulfide_bonds.append('%s %s%s%s' % (line[11:14], line[15:16],
                                                   line[17:21], line[21:22]))
            disulfide_bonds.append('%s %s%s%s' % (line[25:28], line[29:30],
                                                   line[31:35], line[35:36]))

        # Extracts unit cell parameters
        elif line[0:6].strip() == 'CRYST1':
            a = float(line[6:15])
            b = float(line[15:24])
            c = float(line[24:33])
            alpha = math.radians(float(line[33:40]))
            beta = math.radians(float(line[40:47]))
            gamma = math.radians(float(line[47:54]))
            unit_cell_params.extend((a, b, c, alpha, beta, gamma))

        # Extracts non-hydrogen, non-0 occupancy ATOM / HETATM records
        elif (
                line[0:6].strip() in ['ATOM', 'HETATM']
            and line[76:78].strip() != 'H'
            and float(line[54:60].strip()) > 0
            ):
            header = False
            filtered_pdb_lines.append(line)

            # Checks that all disulfide bonds have been refined with 100% occupancy.
            for bond in disulfide_bonds:
                if bond in line and float(line[54:60]) != 1.0:
                    exit = True
                    print('\n\nERROR: One or more disulfide bonds has been '
                          'refined with an occupancy of less than 1.\nTo '
                          'enable damage detection, disulfide bonds should be '
                          'refined as single occupancy\nrather than in '
                          'alternate oxidised and reduced conformations.\n'
                          'Terminating RABDAM run.\n')

            # Checks that all macromolecular atoms in single conformers have an
            # occupancy of 1, and that the occupancies of counterpart atoms in
            # alternate conformers sum to 1.
            if (
                    '{}_{}'.format(line[21:22], line[17:20]) in seqres
                and float(line[54:60]) != 1.0
            ):
                atom_id = '{}_{}_{}_{}'.format(
                    line[21:22], line[22:27].replace(' ', ''),
                    line[17:20].strip(), line[12:16].strip()
                )
                if not atom_id in atom_ids:
                    atom_ids.append(atom_id)
                index = atom_ids.index(atom_id)
                if not line[16:17] in conformers[index]:
                    conformers[index].append(line[16:17])
                    occupancies[index].append(float(line[54:60]))

        elif line[0:6].strip() in ['CONECT', 'MASTER', 'END']:
            footer = True

        if header is True:
            header_lines.append(line)
        elif footer is True:
            footer_lines.append(line)

    orig_pdb.close()

    clean_au_file = ''
    clean_au_list = []
    if exit is False:
        # Completes check for alternate conformer occupancies summing to 1
        for index, atom_id in enumerate(atom_ids):
            if np.sum(occupancies[index]) != 1.0:
                pause = True
                print('ERROR: Atom {} has been refined with sub-1 '
                      'occupancy.'.format(atom_id))
        if pause is True:
            print('')

        # Retains only the most probable alternate conformers. In the case
        # where more than one conformer is equally most probable, (only) that
        # which is listed first in the input PDB file is retained.
        alternate_conformers_chainresnum = []
        alternate_conformers_label = []
        alternate_conformers_occupancy = []
        for line in filtered_pdb_lines:
            if line[16:17].strip() != '':
                alternate_conformers_chainresnum.append(line[21:27].replace(' ', ''))
                alternate_conformers_label.append(line[16:17].strip())
                alternate_conformers_occupancy.append(float(line[54:60].strip()))

        df = pd.DataFrame({'chainresnum': alternate_conformers_chainresnum,
                           'conformer': alternate_conformers_label,
                           'occupancy': alternate_conformers_occupancy})
        df = df.drop_duplicates()
        chainresnum = df['chainresnum'].tolist()
        conformer = df['conformer'].tolist()
        occupancy = df['occupancy'].tolist()

        alternate_conformers = {}
        chainresnum_set = set(chainresnum)
        for number in chainresnum_set:
            indices = []
            a = 'A1'
            b = 'A'
            c = 0
            for index_1, value_1 in enumerate(chainresnum):
                if value_1 == number:
                    indices.append(index_1)
            for index_2 in indices:
                if occupancy[index_2] > c:
                    a = chainresnum[index_2]
                    b = conformer[index_2]
                    c = occupancy[index_2]
            alternate_conformers[a] = b

        for line in filtered_pdb_lines:
            if ((line[16:17].strip() == '')
                or (line[16:17].strip() != ''
                    and alternate_conformers[line[21:27].replace(' ', '')] == line[16:17].strip())
                ):
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
                clean_au_list.append(new_atom)

        clean_au_file = '%s_asymmetric_unit.pdb' % pdb_file_path
        makePDB(header_lines, clean_au_list, footer_lines, clean_au_file, 'Bfactor')

    return (exit, pause, clean_au_file, clean_au_list, header_lines,
            footer_lines, unit_cell_params)


def genPDBCURinputs(PDBCURinputFile):
    # Creates input file for the CCP4 suite program PDBCUR, instructing it to
    # generate the unit cell from an input PDB file of an asymmetric unit

    print('Writing input file for PDBCUR at %s' % PDBCURinputFile)
    input_file = open(PDBCURinputFile, 'w')
    input_file.write('genunit\n')
    input_file.close()


def runPDBCUR(clean_au_file, PDBCURoutputPDB, PDBCURinputFile, PDBCURlog):
    # Runs PDBCUR from the command line with the operations specified in the
    # PDBCUR input file.

    import os

    runPDBCURcommand = 'pdbcur xyzin %s xyzout %s < %s > %s' % (
        clean_au_file, PDBCURoutputPDB, PDBCURinputFile, PDBCURlog
        )
    print('Running PDBCUR (Winn et al. 2011) to generate unit cell')
    os.system(runPDBCURcommand)
    print('PDBCUR log is printed below\n')
    PDBCURlogText = open(PDBCURlog, 'r')
    for line in PDBCURlogText:
        print('%s' % line)
    PDBCURlogText.close()
    os.remove(PDBCURlog)
    os.remove(PDBCURinputFile)

    print('\n################################################################\n'
          '\n################################################################\n'
          '\n################################################################\n')
