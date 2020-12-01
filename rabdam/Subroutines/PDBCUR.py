
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


def make_cryst1_line_from_mmcif(space_group_rec, exit):
    """
    Constructs CRYST1 line in PDB file format (required to run PDBCUR)
    """

    a = None
    b = None
    c = None
    alpha = None
    beta = None
    gamma = None
    sGroup = None
    z = None

    for line in space_group_rec:
        if line.startswith('_cell.length_a '):
            try:
                a = '%9.3f' % float(line.split()[-1])
                if len(a) != 9:
                    exit = True
            except ValueError:
                exit = True
        elif line.startswith('_cell.length_b '):
            try:
                b = '%9.3f' % float(line.split()[-1])
                if len(b) != 9:
                    exit = True
            except ValueError:
                exit = True
        elif line.startswith('_cell.length_c '):
            try:
                c = '%9.3f' % float(line.split()[-1])
                if len(c) != 9:
                    exit = True
            except ValueError:
                exit = True
        elif line.startswith('_cell.angle_alpha '):
            try:
                alpha = '%7.2f' % float(line.split()[-1])
                if len(alpha) != 7:
                    exit = True
            except ValueError:
                exit = True
        elif line.startswith('_cell.angle_beta '):
            try:
                beta = '%7.2f' % float(line.split()[-1])
                if len(beta) != 7:
                    exit = True
            except ValueError:
                exit = True
        elif line.startswith('_cell.angle_gamma '):
            try:
                gamma = '%7.2f' % float(line.split()[-1])
                if len(gamma) != 7:
                    exit = True
            except ValueError:
                exit = True
        elif line.startswith('_symmetry.space_group_name_H-M '):
            sGroup = line.replace('_symmetry.space_group_name_H-M', '')
            sGroup = sGroup.replace("'", '')
            sGroup = sGroup.strip()
            if len(sGroup) > 11:
                exit = True
        elif line.startswith('_cell.Z_PDB '):
            z = line.split()[-1].strip()
            if len(z) > 4:
                exit = True

    if any(x is None for x in [a, b, c, alpha, beta, gamma, sGroup, z]):
        exit = True

    if exit is True:
        print('\n\nERROR: Failed to construct CRYST1 line from input mmCIF file')
        cryst1_line = ''
    else:
        cryst1_line = ''.join(['CRYST1', a, b, c, alpha, beta, gamma, ' ',
                               sGroup.ljust(11), z.rjust(4), '          ', '\n'])

    return cryst1_line, exit


def find_disulfides_from_mmcif(disulfide_rec, exit):
    """
    Extracts info about Cys residues involved in disulfide bonds from input
    mmCIF format file
    """

    disulfide_bonds = {}
    prop_indices = {}
    prop_num = 0
    for line in disulfide_rec:
        if line.startswith('_struct_conn.'):
            prop = line.split('.')[1].strip()
            prop_indices[prop] = prop_num
            prop_num += 1

        elif line.startswith('disulf'):
            try:
                line = line.split()
                disulf_num = int(line[prop_indices['id']].replace('disulf', ''))

                chain1 = line[prop_indices['ptnr1_label_asym_id']]
                resnum1 = int(line[prop_indices['ptnr1_label_seq_id']])
                inscode1 = line[prop_indices['pdbx_ptnr1_PDB_ins_code']]
                res1 = [chain1, resnum1, inscode1]

                chain2 = line[prop_indices['ptnr2_label_asym_id']]
                resnum2 = int(line[prop_indices['ptnr2_label_seq_id']])
                inscode2 = line[prop_indices['pdbx_ptnr2_PDB_ins_code']]
                res2 = [chain2, resnum2, inscode2]

                disulfide_bonds[disulf_num] = [res1, res2]
            except (KeyError, IndexError, ValueError):
                print('Disulfide bond records failed to be parsed by RABDAM - '
                      'check their formatting')
                disulfide_bonds = {}
                exit = True
                break

    return disulfide_bonds, exit


def find_disulfides_from_pdb(disulfide_rec, exit):
    """
    Extracts info about Cys residues involved in disulfide bonds from input
    PDB format file
    """

    disulfide_bonds = {}
    for line in disulfide_rec:
        try:
            disulf_num = int(line[7:10])

            chain1 = line[15:16].strip()
            resnum1 = int(line[17:21])
            inscode1 = line[21:22].strip()
            res1 = [chain1, resnum1, inscode1]

            chain2 = line[29:30].strip()
            resnum2 = int(line[31:35])
            inscode2 = line[35:36].strip()
            res2 = [chain2, resnum2, inscode2]

            disulfide_bonds[disulf_num] = [res1, res2]
        except (IndexError, ValueError):
            print('Disulfide (SSBOND) bond records failed to be parsed by '
                  'RABDAM - check their formatting')
            disulfide_bonds = {}
            exit = True
            break

    return disulfide_bonds, exit


def parse_seqres_from_mmcif(seqres_rec, exit):
    """
    Extracts a list of unique residue identities in the asymmetric unit from an
    input mmCIF format file
    """

    seqres = []
    prop_indices = {}
    prop_num = 0
    for line in seqres_rec:
        if line.startswith('_entity_poly_seq.'):
            prop = line.split('.')[1].strip()
            prop_indices[prop] = prop_num
            prop_num += 1
        elif not any(line.startswith(x) for x in ['loop_', '_entity_poly_seq.']):
            line = line.split()
            try:
                if not line[prop_indices['mon_id']] in seqres:
                    seqres.append(line[prop_indices['mon_id']])
            except (IndexError, KeyError):
                print('Sequence (_entity_poly_seq) records failed to be parsed '
                      'by RABDAM - check their formatting')
                seqres = []
                break

    if seqres == []:
        exit = True

    return seqres, exit


def parse_seqres_from_pdb(seqres_rec, exit):
    """
    Extracts a list of unique residue identities in the asymmetric unit from an
    input PDB format file
    """

    seqres = []
    for line in seqres_rec:
        if len(line) < 19:
            seqres = []
            break
        else:
            res_list = line[19:70].split()
            for res in res_list:
                if not res in seqres:
                    seqres.append(res)

    if seqres == []:
        print('SEQRES records failed to be parsed by RABDAM - check their '
              'formatting')
        exit = True

    return seqres, exit


def parse_atom_rec_from_mmcif(atom_rec, exit):
    """
    Parses ATOM / HETATM records in mmCIF format input file
    """

    import copy

    if __name__ == 'Subroutines.PDBCUR':
        from Subroutines.parsePDB import atom
    else:
        from rabdam.Subroutines.parsePDB import atom

    atoms_list = []
    prop_indices = {}
    prop_num = 0
    for line in atom_rec:
        if line.startswith('_atom_site.'):
            prop = line.split('.')[1].strip()
            prop_indices[prop] = prop_num
            prop_num += 1
        elif line[0:6].strip() in ['ATOM', 'HETATM']:
            try:
                line = line.split()
                new_atom = atom()

                new_atom.lineID = line[prop_indices['group_PDB']]
                new_atom.conformer = line[prop_indices['label_alt_id']].replace('?', '').replace('.', '')
                new_atom.insCode = line[prop_indices['pdbx_PDB_ins_code']].replace('?', '').replace('.', '')
                new_atom.element = line[prop_indices['type_symbol']].replace('?', '').replace('.', '')
                new_atom.charge = line[prop_indices['pdbx_formal_charge']].replace('?', '').replace('.', '')
                new_atom.chainID = line[prop_indices['auth_asym_id']].replace('?', '').replace('.', '')
                new_atom.origChainID = line[prop_indices['label_asym_id']].replace('?', '').replace('.', '')
                new_atom.resiType = line[prop_indices['auth_comp_id']].replace('?', '').replace('.', '')
                new_atom.origResiType = line[prop_indices['label_comp_id']].replace('?', '').replace('.', '')
                new_atom.atomType = line[prop_indices['auth_atom_id']].replace('?', '').replace('.', '')
                if len(new_atom.atomType) > 3 and '\"' in new_atom.atomType:
                    new_atom.atomType = new_atom.atomType.replace('\"', '')
                new_atom.origAtomType = line[prop_indices['label_atom_id']].replace('?', '').replace('.', '')
                new_atom.resiNum = line[prop_indices['auth_seq_id']].replace('?', '').replace('.', '')
                try:
                    new_atom.resiNum = int(new_atom.resiNum)
                except ValueError:
                    print('ERROR: Encountered non-numeric residue number {}'.format(line[prop_indices['auth_seq_id']]))
                    exit = True
                new_atom.origResiNum = line[prop_indices['label_seq_id']].replace('?', '').replace('.', '')
                try:
                    new_atom.atomNum = int(line[prop_indices['id']])
                except ValueError:
                    print('ERROR: Encountered non-numeric atom number {}'.format(line[prop_indices['id']]))
                    exit = True
                try:
                    new_atom.xyzCoords = [[float(line[prop_indices['Cartn_x']])],
                                          [float(line[prop_indices['Cartn_y']])],
                                          [float(line[prop_indices['Cartn_z']])]]
                except ValueError:
                    print('ERROR: Encountered non-numeric coordinate: '
                          '{}'.format([[line[prop_indices['Cartn_x']]],
                                       [line[prop_indices['Cartn_y']]],
                                       [line[prop_indices['Cartn_z']]]]))
                    exit = True
                try:
                    new_atom.occupancy = float(line[prop_indices['occupancy']])
                except ValueError:
                    print('ERROR: Encountered non-numeric occupancy value {}'.format(line[prop_indices['occupancy']]))
                    exit = True
                try:
                    new_atom.bFactor = float(line[prop_indices['B_iso_or_equiv']])
                except ValueError:
                    print('ERROR: Encountered non-numeric B-factor value {}'.format(line[prop_indices['B_iso_or_equiv']]))
                    exit = True

                # Checks all required properties have been assigned to the atom object
                none_props = [
                    new_atom.lineID, new_atom.atomNum, new_atom.atomType,
                    new_atom.conformer, new_atom.resiType, new_atom.chainID,
                    new_atom.resiNum, new_atom.insCode, new_atom.xyzCoords,
                    new_atom.occupancy, new_atom.bFactor, new_atom.element,
                    new_atom.charge, new_atom.origResiNum, new_atom.origResiType,
                    new_atom.origChainID, new_atom.origAtomType
                ]
                if any(prop is None for prop in none_props):
                    exit = True
                else:
                    atoms_list.append(new_atom)

                try:
                    model_num = float(line[prop_indices['pdbx_PDB_model_num']])
                    if not model_num == 1:
                        exit = True
                        print('\n\nERROR: More than one model in the input mmCIF '
                              'file.\nPlease submit a PDB file containing a single '
                              'model for BDamage analysis.\nTerminating RABDAM run.\n')
                except ValueError:
                    print('ERROR: Model number not recognised '
                          '{}'.format(line[prop_indices['pdbx_PDB_model_num']]))
                    exit = True

            except (KeyError, IndexError):
                exit = True
                print('\n\nERROR: Failed to parse ATOM/HETATM records - check '
                      'their formatting\n')

    if exit is True:
        atoms_list = []

    return atoms_list, exit


def parse_atom_rec_from_pdb(atom_rec, exit):
    """
    Parses ATOM / HETATM records in PDB format input file
    """

    if __name__ == 'Subroutines.PDBCUR':
        from Subroutines.parsePDB import atom
    else:
        from rabdam.Subroutines.parsePDB import atom

    atoms_list = []

    for line in atom_rec:
        new_atom = atom()

        new_atom.lineID = line[0:6].strip()
        new_atom.atomType = line[12:16].strip()
        new_atom.origAtomType = line[12:16].strip()
        new_atom.conformer = line[16:17].strip()
        new_atom.resiType = line[17:20].strip()
        new_atom.origResiType = line[17:20].strip()
        new_atom.chainID = line[21:22].strip()
        new_atom.origChainID = line[21:22].strip()
        new_atom.insCode = line[26:27].strip()
        new_atom.element = line[76:78].strip()
        new_atom.charge = line[78:80].strip()

        try:
            new_atom.atomNum = int(line[6:11].strip())
        except ValueError:
            print('ERROR: Encountered non-numeric atom number '
                  '{}'.format(line[6:11].strip()))
        try:
            new_atom.resiNum = int(line[22:26].strip())
            new_atom.origResiNum = line[22:26].strip()
        except ValueError:
            print('ERROR: Encountered non-numeric residue number '
                  '{}'.format(line[22:26].strip()))
        try:
            new_atom.xyzCoords = [[float(line[30:38].strip())],
                                  [float(line[38:46].strip())],
                                  [float(line[46:54].strip())]]
        except ValueError:
            print('ERROR: Encountered non-numeric coordinates: '
                  '{}'.format([[line[30:38].strip()],
                               [line[38:46].strip()],
                               [line[46:54].strip()]]))
        try:
            new_atom.occupancy = float(line[54:60].strip())
        except ValueError:
            print('ERROR: Encountered non-numeric occupancy value '
                  '{}'.format(line[54:60].strip()))
        try:
            new_atom.bFactor = float(line[60:66].strip())
        except ValueError:
            print('ERROR: Encountered non-numeric B-factor value '
                  '{}'.format(line[60:66].strip()))

        # Checks all required properties have been assigned to the atom object
        none_props = [
            new_atom.lineID, new_atom.atomNum, new_atom.atomType,
            new_atom.conformer, new_atom.resiType, new_atom.chainID,
            new_atom.resiNum, new_atom.insCode, new_atom.xyzCoords,
            new_atom.occupancy, new_atom.bFactor, new_atom.element,
            new_atom.charge, new_atom.origResiNum, new_atom.origResiType,
            new_atom.origChainID, new_atom.origAtomType
        ]
        if any(prop is None for prop in none_props):
            exit = True
        else:
            atoms_list.append(new_atom)

    if exit is True:
        atoms_list = []

    return atoms_list, exit


def parse_mmcif_file(pathToInput):
    """
    Parses input mmCIF file
    """

    if __name__ == 'Subroutines.PDBCUR':
        from Subroutines.parsePDB import atom
    else:
        from rabdam.Subroutines.parsePDB import atom

    exit = False

    with open('%s' % pathToInput, 'r') as f:
        input_cif = [subsection for subsection in f.read().split('#')
                     if subsection.strip() != '']

    atom_rec = []
    disulfide_rec = []
    space_group = []
    seqres_rec = []

    for subsection in input_cif:
        if '_atom_site.group_PDB' in subsection:
            atom_rec += [line.strip('\r') for line in subsection.split('\n')
                         if line.strip() != '']
        elif 'disulf1' in subsection:
            disulfide_rec += [line.strip('\r') for line in subsection.split('\n')
                              if not line.strip() in ['', 'loop_']]
        elif ('_cell.length_a' in subsection
             or '_symmetry.space_group_name_H-M' in subsection
        ):
            space_group += [line.strip('\r') for line in subsection.split('\n')
                            if not line.strip() in ['', 'loop_']]
        elif '_entity_poly_seq.' in subsection:
            seqres_rec += [line.strip('\r') for line in subsection.split('\n')
                           if not line.strip() in ['', 'loop_']]

    # Constructs CRYST1 line of PDB file
    cryst1_line, exit = make_cryst1_line_from_mmcif(space_group, exit)
    # Extracts disulfide bonds
    disulfide_bonds, exit = find_disulfides_from_mmcif(disulfide_rec, exit)
    # Extracts SEQRES info from PDB file
    seqres, exit = parse_seqres_from_mmcif(seqres_rec, exit)
    # Extracts ATOM / HETATM lines
    atoms_list, exit = parse_atom_rec_from_mmcif(atom_rec, exit)

    return atoms_list, disulfide_bonds, seqres, cryst1_line, exit


def parse_pdb_file(pathToInput):
    """
    """

    exit = False

    atom_rec = []
    disulfide_rec = []
    seqres_rec = []
    cryst1_line = ''

    with open('%s' % pathToInput, 'r') as f:
        input_pdb = [line.strip('\r') for line in f.read().split('\n')
                     if line.strip() != '']

    for line in input_pdb:
        if len(line) != 80:
            print('\n\nERROR: PDB file formatting incorrect - encountered line '
                  'of less than 80 characters.\nTerminating RABDAM run.\n')
            exit = True
            break
        else:
            if line[0:6].strip() in ['ATOM', 'HETATM']:
                atom_rec.append(line)
            elif line[0:6] == 'SSBOND':
                disulfide_rec.append(line)
            elif line[0:6] == 'SEQRES':
                seqres_rec.append(line)
            elif line[0:6] == 'CRYST1':
                cryst1_line = line
            elif line[0:5] == 'MODEL':
                exit = True
                print('\n\nERROR: More than one model present in input PDB file.\n'
                      'Please submit a PDB file containing a single model for '
                      'BDamage analysis.\n'
                      'Terminating RABDAM run.\n')

    # Extracts disulfide bonds
    disulfide_bonds, exit = find_disulfides_from_pdb(disulfide_rec, exit)
    # Extracts SEQRES info from PDB file
    seqres, exit = parse_seqres_from_pdb(seqres_rec, exit)
    # Extracts ATOM / HETATM lines
    atoms_list, exit = parse_atom_rec_from_pdb(atom_rec, exit)

    return atoms_list, disulfide_bonds, seqres, cryst1_line, exit


def clean_atom_rec(atoms_list, disulfide_bonds, seqres, cryst1_line,
                   file_name_start):
    """
    Filters the ATOM / HETATM records to remove hydrogen atoms and 0 occupancy
    atoms, as well as retaining only the most probable alternate conformers (in
    the case where multiple alternate conformers of a particular atom have the
    same occupancy, only that which is first listed in the PDB file is
    retained). Note that the most probable conformer is selected on a per-atom
    rather than a per-residue basis.
    """

    import copy
    import math
    import shutil
    import numpy as np
    import pandas as pd

    if __name__ == 'Subroutines.PDBCUR':
        from Subroutines.makeDataFrame import makePDB
    else:
        from rabdam.Subroutines.makeDataFrame import makePDB

    # Checks that any disulphide bonds have been refined with 100% occupancy.
    # Then removes anisotropic Bfactors, hydrogen atoms and 0 occupancy atoms.
    exit = False
    pause = False

    filtered_atoms_list = []
    unit_cell_params = []

    atom_ids = {}

    for atm in atoms_list:
        # Extracts non-hydrogen, non-0 B-factor and non-0 occupancy ATOM /
        # HETATM records
        if atm.element != 'H' and atm.bFactor > 0 and 0 < atm.occupancy <= 1:
            filtered_atoms_list.append(copy.deepcopy(atm))

            # Checks that all disulfide bonds between CYS residue pairs have
            # been refined with 100% occupancy.
            for bond_num, res_list in disulfide_bonds.items():
                res1 = res_list[0]
                chain1 = res1[0]
                resnum1 = res1[1]
                inscode1 = res1[2]

                res2 = res_list[1]
                chain2 = res2[0]
                resnum2 = res2[1]
                inscode2 = res2[2]

                if (
                    atm.occupancy != 1
                    and
                    (
                        (atm.chainID == chain1 and atm.resiNum == resnum1
                         and atm.insCode == inscode1 and atm.resiType == 'CYS')
                        or
                        (atm.chainID == chain2 and atm.resiNum == resnum2
                         and atm.insCode == inscode2 and atom.resiType == 'CYS'))
                ):
                    pause = True
                    print('\n\nERROR: One or more disulfide bonds has been '
                          'refined with an occupancy of less than 1.\nTo '
                          'enable damage detection, disulfide bonds should be '
                          'refined as single occupancy\nrather than in '
                          'alternate oxidised and reduced conformations.\n')

            # Checks that all macromolecular atoms in single conformers have an
            # occupancy of 1, and that the occupancies of counterpart atoms in
            # alternate conformers sum to 1.
            if atm.occupancy != 1:
                atom_id = '_'.join([atm.chainID, str(atm.resiNum), atm.insCode, atm.atomType])
                if not atom_id in atom_ids:
                    atom_ids[atom_id] = {}
                atom_ids[atom_id][atm.conformer] = atm.occupancy

    # Completes check for alternate conformer occupancies summing to 1
    discarded_atoms_list = []
    for atom_id, occupancies in atom_ids.items():
        if sum(occupancies.values()) != 1.0:
            pause = True
            print('ERROR: Atom {} has been refined with sub-1 '
                  'occupancy.'.format(atom_id))

        # Single conformer selected on a per-atom basis - hence selected atoms
        # could be a mix of residue conformer "A" and residue conformer "B"
        conformer_index = np.argmax(list(occupancies.values()))
        conformers = list(occupancies.keys())
        discarded_atoms_list += [
            '{}_{}'.format(atom_id, conformer) for index, conformer
            in enumerate(conformers) if not index == conformer_index
        ]

    filtered_atoms_list = [
        atm for atm in filtered_atoms_list if not '_'.join([atm.chainID,
        str(atm.resiNum), atm.insCode, atm.atomType, atm.conformer])
        in discarded_atoms_list
    ]
    if len(filtered_atoms_list) == 0:
        exit = True
        print('\n\nERROR: No atoms retained for BDamage analysis after cleaning'
              ' input file')

    clean_au_file = '%s_asymmetric_unit.pdb' % file_name_start
    pdb_fail = makePDB([cryst1_line], filtered_atoms_list, [], seqres, clean_au_file, 'bfactor')
    if pdb_fail is True:
        exit  = True
        print('\n\nERROR: Failed to make input PDB file to feed into PDBCUR.\n'
              'Check formatting of input PDB/mmCIF file input file')

    return exit, pause, filtered_atoms_list, clean_au_file


def genPDBCURinputs(PDBCURinputFile):
    """
    Creates input file for the CCP4 suite program PDBCUR, instructing it to
    generate the unit cell from an input PDB file of an asymmetric unit
    """

    print('Writing input file for PDBCUR at %s' % PDBCURinputFile)
    input_file = open(PDBCURinputFile, 'w')
    input_file.write('genunit\n')
    input_file.close()


def runPDBCUR(clean_au_file, PDBCURoutputPDB, PDBCURinputFile, PDBCURlog):
    """
    Runs PDBCUR from the command line with the operations specified in the
    PDBCUR input file.
    """

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
