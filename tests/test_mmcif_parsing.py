
# RABDAM
# Copyright (C) 2024 Garman Group, University of Oxford

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


# python -m unittest tests/test_mmcif_parsing.py

import unittest
from rabdam.Subroutines.parsePDB import atom, download_mmcif
from rabdam.Subroutines.PDBCUR import (
    get_res_rfree_temp_from_mmcif, make_cryst1_line_from_mmcif,
    parse_atom_rec_from_mmcif
)
from rabdam.Subroutines.output import write_output_cif, write_all_carbon_cif
from tests.make_atom_rec import gen_atom_objs_list, make_exp_mmcif_rec


class TestClass(unittest.TestCase):

    def test_download_mmcif(self):
        """
        Checks that downloads mmCIF files correctly, and raises the appropriate
        errors e.g. if the file does not exist
        """

        import os
        import shutil

        if not os.path.isdir('tests/temp_files'):
            os.mkdir('tests/temp_files')

        exp_results = {'2BN1': ['tests/temp_files', False],
                       '2BN3': ['tests/tmp_files', True],
                       '1ABC': ['tests/temp_files', True]}
        for pdb, exp_vals in exp_results.items():
            out_dir = exp_vals[0]
            exp_exit = exp_vals[1]
            act_exit = download_mmcif(pdb, '{}/{}.cif'.format(out_dir, pdb))
            self.assertEqual(exp_exit, act_exit)
            if pdb == '2BN1':
                self.assertTrue(os.path.isfile('tests/temp_files/2BN1.cif'))

        shutil.rmtree('tests/temp_files/')

    def test_make_cryst1_line_from_mmcif(self):
        """
        Checks that CRYST1 line (found in PDB format files, and required for
        translating the unit cell into a 3x3 parallelepiped) is generated
        correctly from an mmCIF format input file
        """

        import copy

        # Taken from 2BN1
        space_group_info_1 = [
            '_cell.entry_id           2BN1',
            '_cell.length_a           77.900',
            '_cell.length_b           77.900',
            '_cell.length_c           77.900',
            '_cell.angle_alpha        90.00',
            '_cell.angle_beta         90.00',
            '_cell.angle_gamma        90.00',
            '_cell.Z_PDB              24',
            '_cell.pdbx_unique_axis   ?',
            '_symmetry.entry_id                         2BN1',
            '_symmetry.space_group_name_H-M             \'I 21 3\'',
            '_symmetry.pdbx_full_space_group_name_H-M   ?',
            '_symmetry.cell_setting                     ?',
            '_symmetry.Int_Tables_number                199'
        ]

        # Check CRYST1 line is constructed correctly
        exp_cryst1_line_1 = 'CRYST1   77.900   77.900   77.900  90.00  90.00  90.00 I 21 3                   \n'
        act_cryst1_line_1, exit_1 = make_cryst1_line_from_mmcif(space_group_info_1, False)
        self.assertEqual(exp_cryst1_line_1, act_cryst1_line_1)
        self.assertFalse(exit_1)  

        # Check raises error if dimension is too long
        space_group_info_2 = copy.deepcopy(space_group_info_1)
        space_group_info_2[1] = '_cell.length_a       123477.900'
        exp_cryst1_line_2 = ''
        act_cryst1_line_2, exit_2 = make_cryst1_line_from_mmcif(space_group_info_2, False)
        self.assertEqual(exp_cryst1_line_2, act_cryst1_line_2)
        self.assertTrue(exit_2)

        # Check raises error if dimension is not a numeric value
        space_group_info_3 = copy.deepcopy(space_group_info_1)
        space_group_info_3[2] = '_cell.length_b                A'
        exp_cryst1_line_3 = ''
        act_cryst1_line_3, exit_3 = make_cryst1_line_from_mmcif(space_group_info_3, False)
        self.assertEqual(exp_cryst1_line_3, act_cryst1_line_3)
        self.assertTrue(exit_3)

        # Check CRYST1 line is constructed correctly even though spacing is
        # non-standard. Also check dimension is rounded to the correct number of
        # decimal places (3) in the CRYST1 line.
        space_group_info_4 = copy.deepcopy(space_group_info_1)
        space_group_info_4[3] = '_cell.length_c 18.5'
        exp_cryst1_line_4 = 'CRYST1   77.900   77.900   18.500  90.00  90.00  90.00 I 21 3                   \n'
        act_cryst1_line_4, exit_4 = make_cryst1_line_from_mmcif(space_group_info_4, False)
        self.assertEqual(exp_cryst1_line_4, act_cryst1_line_4)
        self.assertFalse(exit_4)

        # Check CRYST1 line is constructed correctly even though spacing is
        # non-standard.
        space_group_info_5 = copy.deepcopy(space_group_info_1)
        space_group_info_5[4] = '_cell.angle_alpha               15.00'
        exp_cryst1_line_5 = 'CRYST1   77.900   77.900   77.900  15.00  90.00  90.00 I 21 3                   \n'
        act_cryst1_line_5, exit_5 = make_cryst1_line_from_mmcif(space_group_info_5, False)
        self.assertEqual(exp_cryst1_line_5, act_cryst1_line_5)
        self.assertFalse(exit_5)

        # Check CRYST1 line is constructed correctly, and to round angle to
        # the expected number of decimal places (2).
        space_group_info_6 = copy.deepcopy(space_group_info_1)
        space_group_info_6[5] = '_cell.angle_beta         15.5678'
        exp_cryst1_line_6 = 'CRYST1   77.900   77.900   77.900  90.00  15.57  90.00 I 21 3                   \n'
        act_cryst1_line_6, exit_6 = make_cryst1_line_from_mmcif(space_group_info_6, False)
        self.assertEqual(exp_cryst1_line_6, act_cryst1_line_6)
        self.assertFalse(exit_6)

        # Check raises error if angle is too large (should be 7 characters, with
        # 2 characters reserved for after the decimal point)
        space_group_info_7 = copy.deepcopy(space_group_info_1)
        space_group_info_7[5] = '_cell.angle_gamma         123456'
        exp_cryst1_line_7 = ''
        act_cryst1_line_7, exit_7 = make_cryst1_line_from_mmcif(space_group_info_7, False)
        self.assertEqual(exp_cryst1_line_7, act_cryst1_line_7)
        self.assertTrue(exit_7)

        # Check raises error because H-M space group name is longer than 11
        # characters
        space_group_info_8 = copy.deepcopy(space_group_info_1)
        space_group_info_8[7] = '_symmetry.space_group_name_H-M \'123456789101\''
        exp_cryst1_line_8 = ''
        act_cryst1_line_8, exit_8 = make_cryst1_line_from_mmcif(space_group_info_8, False)
        self.assertEqual(exp_cryst1_line_8, act_cryst1_line_8)
        self.assertTrue(exit_8)

        # Check raises error because Z_PDB (= the number of polymeric chains in
        # a unit cell) is longer than 4 characters
        space_group_info_9 = copy.deepcopy(space_group_info_1)
        space_group_info_9[8] = '_cell.Z_PDB 12345'
        exp_cryst1_line_9 = 'CRYST1   77.900   77.900   77.900  90.00  90.00  90.00 I 21 3                   \n'
        act_cryst1_line_9, exit_9 = make_cryst1_line_from_mmcif(space_group_info_9, False)
        self.assertEqual(exp_cryst1_line_9, act_cryst1_line_9)
        self.assertFalse(exit_9)

        # Check raises error because one of the expected unit cell dimensions is
        # not supplied
        space_group_info_10 = copy.deepcopy(space_group_info_1)
        space_group_info_10[1] = ''
        exp_cryst1_line_10 = ''
        act_cryst1_line_10, exit_10 = make_cryst1_line_from_mmcif(space_group_info_10, False)
        self.assertEqual(exp_cryst1_line_10, act_cryst1_line_10)
        self.assertTrue(exit_10)

    def test_get_res_rfree_temp_from_mmcif(self):
        """
        Checks that resolution, Rfree and temperature values are parsed as
        expected from an input mmCIF file
        """

        import copy

        remark_dict_1 = {
            'resolution': '_refine.ls_d_res_high                            1.40 \n_refine.ls_percent_reflns_obs                    99.8 \n',
            'rfree': '_refine.ls_R_factor_R_free                       0.165 \n_refine.ls_R_factor_R_free_error                 ? \n',
            'temperature': '_diffrn.ambient_temp           100 \n_diffrn.ambient_temp_details   ? \n'
        }

        # Check that mmCIF file with no errors (from 2BN1) is parsed correctly
        act_resolution_1, act_rfree_1, act_temp_1 = get_res_rfree_temp_from_mmcif(remark_dict_1)
        exp_resolution_1 = 1.4
        exp_rfree_1 = 0.165
        exp_temp_1 = 100.0
        self.assertEqual(act_resolution_1, exp_resolution_1)
        self.assertEqual(act_rfree_1, exp_rfree_1)
        self.assertEqual(act_temp_1, exp_temp_1)

        # Check that if resolution is a non-numeric value, resolution is set to
        # None
        remark_dict_2 = copy.deepcopy(remark_dict_1)
        remark_dict_2['resolution'] = '1.4.1'
        act_resolution_2, act_rfree_2, act_temp_2 = get_res_rfree_temp_from_mmcif(remark_dict_2)
        exp_resolution_2 = None
        exp_rfree_2 = 0.165
        exp_temp_2 = 100.0
        self.assertEqual(act_resolution_2, exp_resolution_2)
        self.assertEqual(act_rfree_2, exp_rfree_2)
        self.assertEqual(act_temp_2, exp_temp_2)

        # Check that if resolution is a non-numeric value, resolution is set to
        # None
        remark_dict_3 = copy.deepcopy(remark_dict_1)
        remark_dict_3['rfree'] = ''
        act_resolution_3, act_rfree_3, act_temp_3 = get_res_rfree_temp_from_mmcif(remark_dict_3)
        exp_resolution_3 = 1.4
        exp_rfree_3 = None
        exp_temp_3 = 100.0
        self.assertEqual(act_resolution_3, exp_resolution_3)
        self.assertEqual(act_rfree_3, exp_rfree_3)
        self.assertEqual(act_temp_3, exp_temp_3)

        # Check that if resolution is a non-numeric value, resolution is set to
        # None
        remark_dict_4 = copy.deepcopy(remark_dict_1)
        remark_dict_4['temperature'] = '100K'
        act_resolution_4, act_rfree_4, act_temp_4 = get_res_rfree_temp_from_mmcif(remark_dict_4)
        exp_resolution_4 = 1.4
        exp_rfree_4 = 0.165
        exp_temp_4 = None
        self.assertEqual(act_resolution_4, exp_resolution_4)
        self.assertEqual(act_rfree_4, exp_rfree_4)
        self.assertEqual(act_temp_4, exp_temp_4)


    def test_parse_atomrec_from_mmcif(self):
        """
        Checks that ATOM/HETATM records in mmCIF format files are parsed as
        expected
        """

        import copy
        import os
        import shutil

        # Load example mmCIF records
        exp_atom_rec_obj_1, exp_atom_rec_str_1 = make_exp_mmcif_rec()
        if not os.path.isdir('tests/temp_files/'):
            os.mkdir('tests/temp_files/')
        with open('tests/temp_files/input.cif', 'w') as f:
            f.write('data_1ABC\n#\nloop_\n')
            f.write('\n'.join(exp_atom_rec_str_1))
            f.write('\n#')

        # Check that atom records are parsed correctly
        act_atom_rec_1, exit_1 = parse_atom_rec_from_mmcif(
            exp_atom_rec_str_1, 'tests/temp_files/input.cif', False, False
        )
        for index in range(len(act_atom_rec_1)):
            self.assertEqual(exp_atom_rec_obj_1[index], act_atom_rec_1[index])
        self.assertFalse(exit_1)

        # Check raises error if resiNum is not a numeric value
        exp_atom_rec_str_2 = copy.deepcopy(exp_atom_rec_str_1)
        exp_atom_rec_str_2[21] = (
            exp_atom_rec_str_2[21][:74] + '?' + exp_atom_rec_str_2[21][75:]
        )
        act_atom_rec_2, exit_2 = parse_atom_rec_from_mmcif(
            exp_atom_rec_str_2, 'tests/temp_files/input.cif', False, False
        )
        self.assertTrue(exit_2)

        # Check raises error if atomNum is not a numeric value
        exp_atom_rec_str_3 = copy.deepcopy(exp_atom_rec_str_1)
        exp_atom_rec_str_3[21] = (
            exp_atom_rec_str_3[21][:7] + '?' + exp_atom_rec_str_3[21][8:]
        )
        act_atom_rec_3, exit_3 = parse_atom_rec_from_mmcif(
            exp_atom_rec_str_3, 'tests/temp_files/input.cif', False, False
        )
        self.assertTrue(exit_3)

        # Check raises error if xyzcoords are not numeric values
        exp_atom_rec_str_4 = copy.deepcopy(exp_atom_rec_str_1)
        exp_atom_rec_str_4[21] = (
            exp_atom_rec_str_4[21][:36] + '   ?   ' + exp_atom_rec_str_4[21][43:]
        )
        act_atom_rec_4, exit_4 = parse_atom_rec_from_mmcif(
            exp_atom_rec_str_4, 'tests/temp_files/input.cif', False, False
        )
        self.assertTrue(exit_4)

        # Check raises error if occupancy is not a numeric value
        exp_atom_rec_str_5 = copy.deepcopy(exp_atom_rec_str_1)
        exp_atom_rec_str_5[21] = (
            exp_atom_rec_str_5[21][:60] + ' ?  ' + exp_atom_rec_str_5[21][64:]
        )
        act_atom_rec_5, exit_5 = parse_atom_rec_from_mmcif(
            exp_atom_rec_str_5, 'tests/temp_files/input.cif', False, False
        )
        self.assertTrue(exit_5)

        # Check raises error if B-factor is not a numeric value
        exp_atom_rec_str_6 = copy.deepcopy(exp_atom_rec_str_1)
        exp_atom_rec_str_6[21] = (
            exp_atom_rec_str_6[21][:65] + '  ?  ' + exp_atom_rec_str_6[21][70:]
        )
        act_atom_rec_6, exit_6 = parse_atom_rec_from_mmcif(
            exp_atom_rec_str_6, 'tests/temp_files/input.cif', False, False
        )
        self.assertTrue(exit_6)

        # Check raises error if model number is not a numeric value
        exp_atom_rec_str_7 = copy.deepcopy(exp_atom_rec_str_1)
        exp_atom_rec_str_7[21] = (
            exp_atom_rec_str_7[21][:90] + '?' + exp_atom_rec_str_7[21][91:]
        )
        act_atom_rec_7, exit_7 = parse_atom_rec_from_mmcif(
            exp_atom_rec_str_7, 'tests/temp_files/input.cif', False, False
        )
        self.assertTrue(exit_7)

        # Check raises error if model number is not 1
        exp_atom_rec_str_8 = copy.deepcopy(exp_atom_rec_str_1)
        exp_atom_rec_str_8[21] = (
            exp_atom_rec_str_8[21][:90] + '2' + exp_atom_rec_str_8[21][91:]
        )
        act_atom_rec_8, exit_8 = parse_atom_rec_from_mmcif(
            exp_atom_rec_str_8, 'tests/temp_files/input.cif', False, False
        )
        self.assertTrue(exit_8)

        # Check raises error if atom property remains undefined
        exp_atom_rec_str_9 = copy.deepcopy(exp_atom_rec_str_1)
        exp_atom_rec_str_9 = exp_atom_rec_str_9[:2] + exp_atom_rec_str_9[3:]
        act_atom_rec_9, exit_9 = parse_atom_rec_from_mmcif(
            exp_atom_rec_str_9, 'tests/temp_files/input.cif', False, False
        )
        self.assertTrue(exit_9)

        shutil.rmtree('tests/temp_files/')

    def test_write_mmcif_file(self):
        """
        Tests that atom objects are correctly converted into mmCIF format
        """

        import copy
        import numpy as np
        import os
        import shutil

        if not os.path.isdir('tests/temp_files/'):
            os.mkdir('tests/temp_files/')

        # Generate list of atom objects
        atoms_list = gen_atom_objs_list()

        exp_cif_lines_1 = [
            'data_tests/temp_files/test',
            '#',
            'loop_',
            '_atom_site.group_PDB',
            '_atom_site.id',
            '_atom_site.type_symbol',
            '_atom_site.label_atom_id',
            '_atom_site.label_alt_id',
            '_atom_site.label_comp_id',
            '_atom_site.label_asym_id',
            '_atom_site.label_entity_id',
            '_atom_site.label_seq_id',
            '_atom_site.pdbx_PDB_ins_code',
            '_atom_site.Cartn_x',
            '_atom_site.Cartn_y',
            '_atom_site.Cartn_z',
            '_atom_site.occupancy',
            '_atom_site.B_iso_or_equiv',
            '_atom_site.B_Damage',
            '_atom_site.packing_density',
            '_atom_site.average_B_factor',
            '_atom_site.pdbx_formal_charge',
            '_atom_site.auth_seq_id',
            '_atom_site.auth_comp_id',
            '_atom_site.auth_asym_id',
            '_atom_site.auth_atom_id',
            '_atom_site.pdbx_PDB_model_num',
            'ATOM   2    C CA . GLY A 1 1 ? 13.602 45.768 30.728 1.0 19.08 1.3  19 26.2 ? 1    GLY A CA 1.0',
            'HETATM 438  O O  . HOH C 2 . ? 13.128 46.298 34.342 1.0 33.1  0.97 30 25.6 ? 2001 HOH B O  1.1',
            'HETATM 3546 P PA X FAD B 3 . ? 29.643 29.7   51.402 1.0 20.52 2.0  26 20.5 ? 700  FAD B PA 12 ',
            '#',
            ''
        ]

        # Check mmCIF file is created as expected with BDamage entries
        write_output_cif(atoms_list, 'tests/temp_files/test', True)
        with open('tests/temp_files/test.cif', 'r') as f:
            act_cif_lines_1 = f.read().split('\n')
        self.assertEqual(exp_cif_lines_1, act_cif_lines_1)

        # Check mmCIF file is created as expected without BDamage entries
        write_output_cif(atoms_list, 'tests/temp_files/test', False)
        exp_cif_lines_2 = []
        bdam_header_lines = [
            '_atom_site.B_Damage', '_atom_site.packing_density',
            '_atom_site.average_B_factor'
        ]
        for line in copy.deepcopy(exp_cif_lines_1):
            if line in bdam_header_lines:
                continue
            elif any(line.startswith(start) for start in ['ATOM', 'HETATM']):
                line = line[:61] + line[74:]
                exp_cif_lines_2.append(line)
            else:
                exp_cif_lines_2.append(line)

        with open('tests/temp_files/test.cif', 'r') as f:
            act_cif_lines_2 = f.read().split('\n')
        self.assertEqual(exp_cif_lines_2, act_cif_lines_2)

        # Check all-carbon mmCIF creation
        atom_id_list = [15, 851, 35]
        atoms_list_3 = np.array([
            [13.602, 45.768, 30.728],
            [13.128, 46.298, 34.342],
            [29.643, 29.7, 51.402]
        ])
        write_all_carbon_cif(atoms_list_3, atom_id_list, 'tests/temp_files/test')
        exp_cif_lines_3 = [
            'data_tests/temp_files/test',
            '#',
            'loop_',
            '_atom_site.group_PDB',
            '_atom_site.id',
            '_atom_site.type_symbol',
            '_atom_site.label_atom_id',
            '_atom_site.label_alt_id',
            '_atom_site.label_comp_id',
            '_atom_site.label_asym_id',
            '_atom_site.label_entity_id',
            '_atom_site.label_seq_id',
            '_atom_site.pdbx_PDB_ins_code',
            '_atom_site.Cartn_x',
            '_atom_site.Cartn_y',
            '_atom_site.Cartn_z',
            '_atom_site.occupancy',
            '_atom_site.B_iso_or_equiv',
            '_atom_site.pdbx_formal_charge',
            '_atom_site.auth_seq_id',
            '_atom_site.auth_comp_id',
            '_atom_site.auth_asym_id',
            '_atom_site.auth_atom_id',
            '_atom_site.pdbx_PDB_model_num',
            'ATOM 15  C C ? GLY A 1 1 ? 13.602 45.768 30.728 1.00 50.0 ? 1 GLY A C 1',
            'ATOM 851 C C ? GLY A 1 1 ? 13.128 46.298 34.342 1.00 50.0 ? 1 GLY A C 1',
            'ATOM 35  C C ? GLY A 1 1 ? 29.643 29.7   51.402 1.00 50.0 ? 1 GLY A C 1',
            '#',
            ''
        ]
        with open('tests/temp_files/test.cif', 'r') as f:
            act_cif_lines_3 = f.read().split('\n')
        self.assertEqual(exp_cif_lines_3, act_cif_lines_3)

        shutil.rmtree('tests/temp_files/')
