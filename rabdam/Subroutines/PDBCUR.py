
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


def parse_atom_rec_from_mmcif(atom_rec, input_cif, exit):
    """
    Parses ATOM / HETATM records in mmCIF format input file
    """

    from iotbx.data_manager import DataManager
    import numpy as np

    if __name__ == 'Subroutines.PDBCUR':
        from Subroutines.parsePDB import atom
    else:
        from rabdam.Subroutines.parsePDB import atom

    dm = DataManager()
    all_models = dm.get_model(input_cif)

    model = all_models.get_hierarchy().models()[0]
    model_atoms = list(model.atoms())
    model_atoms_dict = {}
    for atm in model_atoms:
        model_atoms_dict[' '.join([str(n) for n in list(atm.xyz)])] = atm

    res_is_prot_dict = {}
    res_is_na_dict = {}
    chain_len_dict = {}
    atoms_list = []
    prop_indices = {}
    prop_num = 0
    for line in atom_rec:
        if line.startswith('_atom_site.'):
            prop = line.split('.')[1].strip()
            prop_indices[prop] = prop_num
            prop_num += 1
        elif any(x in line for x in ['ATOM', 'HETATM']):
            try:
                line = line.split()
                new_atom = atom()

                new_atom.lineID = line[prop_indices['group_PDB']]
                new_atom.conformer = line[prop_indices['label_alt_id']].replace('?', '').replace('.', '')
                new_atom.insCode = line[prop_indices['pdbx_PDB_ins_code']].replace('?', '').replace('.', '')
                new_atom.element = line[prop_indices['type_symbol']].replace('?', '').replace('.', '')
                new_atom.charge = ''  # Avoids errors when using cctbx to parse clean_au_file to generate a copy of the unit cell
                #new_atom.charge = line[prop_indices['pdbx_formal_charge']].replace('?', '').replace('.', '')
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
                try:
                    atom_coords = ' '.join([str(n) for n in np.array(new_atom.xyzCoords).flatten()])
                    if not atom_coords in res_is_prot_dict.keys():
                        atm = model_atoms_dict[atom_coords]
                        res_is_prot_dict[atom_coords] = atm.parent().parent().conformers()[0].is_protein()
                    new_atom.protein = res_is_prot_dict[atom_coords]
                except (KeyError, AttributeError):
                    print('ERROR: Unable to identify if atom {} is in a protein'
                        ' chain'.format(line[prop_indices['id']]))
                    exit = True
                try:
                    if not atom_coords in res_is_na_dict.keys():
                        atm = model_atoms_dict[atom_coords]
                        res_is_na_dict[atom_coords] = atm.parent().parent().conformers()[0].is_na()
                    new_atom.na = res_is_na_dict[atom_coords]
                except (KeyError, AttributeError):
                    print('ERROR: Unable to identify if atom {} is in a nucleic acid'
                        ' chain'.format(line[prop_indices['id']]))
                    exit = True
                try:
                    if not atom_coords in chain_len_dict.keys():
                        atm = model_atoms_dict[atom_coords]
                        chain_len_dict[atom_coords] = atm.parent().parent().parent().residue_groups_size()
                    new_atom.chain_len = chain_len_dict[atom_coords]
                except (KeyError, AttributeError):
                    print('ERROR: Unable to identify the length of the parent '
                          'chain of atom {}'.format(line[prop_indices['id']]))
                    exit = True

                # Checks all required properties have been assigned to the atom object
                none_props = [
                    new_atom.lineID, new_atom.atomNum, new_atom.atomType,
                    new_atom.conformer, new_atom.resiType, new_atom.chainID,
                    new_atom.resiNum, new_atom.insCode, new_atom.xyzCoords,
                    new_atom.occupancy, new_atom.bFactor, new_atom.element,
                    new_atom.charge, new_atom.origResiNum, new_atom.origResiType,
                    new_atom.origChainID, new_atom.origAtomType,
                    new_atom.protein, new_atom.na, new_atom.chain_len
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


def parse_atom_rec_from_pdb(atom_rec, input_pdb, exit):
    """
    Parses ATOM / HETATM records in PDB format input file
    """

    from iotbx.data_manager import DataManager
    import numpy as np

    if __name__ == 'Subroutines.PDBCUR':
        from Subroutines.parsePDB import atom
    else:
        from rabdam.Subroutines.parsePDB import atom

    dm = DataManager()
    all_models = dm.get_model(input_pdb)

    model = all_models.get_hierarchy().models()[0]
    model_atoms = list(model.atoms())
    model_atoms_dict = {}
    for atm in model_atoms:
        model_atoms_dict[' '.join([str(n) for n in list(atm.xyz)])] = atm

    res_is_prot_dict = {}
    res_is_na_dict = {}
    chain_len_dict = {}
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
        new_atom.charge = ''  # Avoids errors when using cctbx to parse
        # clean_au_file to generate a copy of the unit cell
        """
        try:
            new_atom.charge = line[78:80].strip()
        except IndexError:
            new_atom.charge = ''
        """
        try:
            atom_coords = ' '.join([str(n) for n in np.array(new_atom.xyzCoords).flatten()])
            if not atom_coords in res_is_prot_dict.keys():
                atm = model_atoms_dict[atom_coords]
                res_is_prot_dict[atom_coords] = atm.parent().parent().conformers()[0].is_protein()
            new_atom.protein = res_is_prot_dict[atom_coords]
        except (KeyError, AttributeError):
            print('ERROR: Unable to identify if atom {} is in a protein'
                ' chain'.format(line[6:11].strip()))
            exit = True
        try:
            if not atom_coords in res_is_na_dict.keys():
                atm = model_atoms_dict[atom_coords]
                res_is_na_dict[atom_coords] = atm.parent().parent().conformers()[0].is_na()
            new_atom.na = res_is_na_dict[atom_coords]
        except (KeyError, AttributeError):
            print('ERROR: Unable to identify if atom {} is in a nucleic acid'
                ' chain'.format(line[6:11].strip()))
            exit = True
        try:
            if not atom_coords in chain_len_dict.keys():
                atm = model_atoms_dict[atom_coords]
                chain_len_dict[atom_coords] = atm.parent().parent().parent().residue_groups_size()
            new_atom.chain_len = chain_len_dict[atom_coords]
        except (KeyError, AttributeError):
            print('ERROR: Unable to identify the length of the parent chain of '
                  'atom {}'.format(line[6:11].strip()))
            exit = True

        # Checks all required properties have been assigned to the atom object
        none_props = [
            new_atom.lineID, new_atom.atomNum, new_atom.atomType,
            new_atom.conformer, new_atom.resiType, new_atom.chainID,
            new_atom.resiNum, new_atom.insCode, new_atom.xyzCoords,
            new_atom.occupancy, new_atom.bFactor, new_atom.element,
            new_atom.charge, new_atom.origResiNum, new_atom.origResiType,
            new_atom.origChainID, new_atom.origAtomType, new_atom.protein,
            new_atom.na, new_atom.chain_len
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
    space_group = []

    for subsection in input_cif:
        if '_atom_site.group_PDB' in subsection:
            atom_rec += [line.strip('\r') for line in subsection.split('\n')
                         if line.strip() != '']
        elif ('_cell.length_a' in subsection
             or '_symmetry.space_group_name_H-M' in subsection
        ):
            space_group += [line.strip('\r') for line in subsection.split('\n')
                            if not line.strip() in ['', 'loop_']]

    # Constructs CRYST1 line of PDB file
    cryst1_line, exit = make_cryst1_line_from_mmcif(space_group, exit)
    # Extracts ATOM / HETATM lines
    atoms_list, exit = parse_atom_rec_from_mmcif(atom_rec, pathToInput, exit)

    return atoms_list, cryst1_line, exit


def parse_pdb_file(pathToInput):
    """
    """

    exit = False

    atom_rec = []
    cryst1_line = ''

    with open('%s' % pathToInput, 'r') as f:
        input_pdb = [line.strip('\r') for line in f.read().split('\n')
                     if line.strip() != '']

    for line in input_pdb:
        if line[0:6].strip() in ['ATOM', 'HETATM']:
            atom_rec.append(line)
        elif line[0:6] == 'CRYST1':
            cryst1_line = line
        elif line[0:5] == 'MODEL':
            exit = True
            print('\n\nERROR: More than one model present in input PDB file.\n'
                    'Please submit a PDB file containing a single model.\n'
                    'Terminating RABDAM run.\n')

    # Extracts ATOM / HETATM lines
    atoms_list, exit = parse_atom_rec_from_pdb(atom_rec, pathToInput, exit)

    return atoms_list, cryst1_line, exit


def clean_atom_rec(atoms_list, cryst1_line, file_name_start):
    """
    Filters the ATOM / HETATM records to remove hydrogen atoms and 0 occupancy
    atoms, as well as retaining only the most probable alternate conformers (in
    the case where multiple alternate conformers of a particular atom have the
    same occupancy, only that which is first listed in the PDB file is
    retained). Note that the most probable conformer is selected on a per-atom
    rather than a per-residue basis.
    """

    import copy
    import numpy as np

    if __name__ == 'Subroutines.PDBCUR':
        from Subroutines.makeDataFrame import makePDB
    else:
        from rabdam.Subroutines.makeDataFrame import makePDB

    exit = False
    pause = False

    filtered_atoms_list = []
    atom_ids = {}

    for atm in atoms_list:
        # Extracts non-hydrogen, non-0 B-factor and non-0 occupancy ATOM /
        # HETATM records
        if atm.element != 'H' and atm.bFactor > 0 and 0 < atm.occupancy <= 1:
            filtered_atoms_list.append(copy.deepcopy(atm))

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
    pdb_fail = makePDB([cryst1_line], filtered_atoms_list, [], clean_au_file, 'bfactor')
    if pdb_fail is True:
        exit  = True
        print('\n\nERROR: Failed to make asymmetric unit PDB file.\n'
              'Check formatting of input PDB/mmCIF file input file')

    return exit, pause, filtered_atoms_list, clean_au_file


def gen_unit_cell(clean_au_file):
    """
    Only need xyz coordinates for unit cell atoms => should speed up the code considerably
    """

    import numpy as np
    from iotbx.data_manager import DataManager

    dm = DataManager()
    model = dm.get_model(clean_au_file)
    xray_structure = model.get_xray_structure()
    p1 = xray_structure.expand_to_p1()
    unit_cell_coords = np.array(p1.sites_cart())

    return unit_cell_coords, [n for n in range(1, (unit_cell_coords.shape[0] + 1))]




