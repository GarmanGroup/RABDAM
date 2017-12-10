
# RABDAM
# Copyright (C) 2017 Garman Group, University of Oxford

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


def clean_pdb_file(pathToPDB, PDBdirectory, batchRun, pdb_file_path):
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
    import pandas as pd
    from parsePDB import atom

    # Removes anisotropic Bfactors, hydrogen atoms and 0 occupancy atoms
    multi_model = False
    filtered_pdb_lines = []
    header_lines = []
    footer_lines = []
    unit_cell_params = []
    orig_pdb = open('%s' % pathToPDB, 'r')
    header = True
    footer = False
    for line in orig_pdb:
        if line.replace(' ', '').startswith('MODEL2'):
            multi_model = True
            if batchRun is False:
                shutil.rmtree('%s' % PDBdirectory)
                sys.exit('\n\nMore than one model present in input PDB file.\n'
                         'Please submit a PDB file containing a single model '
                         'for BDamage analysis.\n')
        elif line[0:6].strip() == 'CRYST1':
            params = line.split()
            a = float(params[1])
            b = float(params[2])
            c = float(params[3])
            alpha = math.radians(float(params[4]))
            beta = math.radians(float(params[5]))
            gamma = math.radians(float(params[6]))
            unit_cell_params.extend((a, b, c, alpha, beta, gamma))
        elif (line[0:6].strip() in ['ATOM', 'HETATM']
            and line[76:78].strip() != 'H'
            and float(line[54:60].strip()) > 0
             ):
            filtered_pdb_lines.append(line)
            header = False
        elif line[0:6].strip() == 'CONECT':
            footer = True

        if header is True:
            header_lines.append(line)
        elif footer is True:
            footer_lines.append(line)
    orig_pdb.close()

    # Retains only the most probable alternate conformers. In the case where
    # more than one conformer is equally most probable, (only) that which is
    # listed first in the input PDB file is retained.
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

    clean_au_list = []
    clean_au_file_name = '%s_asymmetric_unit.pdb' % pdb_file_path
    clean_au_file = open('%s' % clean_au_file_name, 'w')

    for line in header_lines:
        clean_au_file.write(line)

    for line in filtered_pdb_lines:
        if ((line[16:17].strip() == '')
            or (line[16:17].strip() != ''
                and alternate_conformers[line[21:27].replace(' ', '')] == line[16:17].strip())
            ):
            clean_au_file.write(line)
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

    for line in footer_lines:
        clean_au_file.write(line)

    clean_au_file.close()

    return (multi_model, clean_au_file_name, clean_au_list, header_lines,
            footer_lines, unit_cell_params)


def genPDBCURinputs(PDBCURinputFile):
    # Creates input file for the CCP4 suite program PDBCUR, instructing it to
    # create the unit cell from an input pdb file of an asymmetric unit

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
