
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


# python -m unittest tests/test_mmcif_parsing.py

import unittest
from rabdam.Subroutines.parsePDB import atom, download_mmcif
from rabdam.Subroutines.PDBCUR import (
    make_cryst1_line_from_mmcif, parse_atom_rec_from_mmcif
)
from rabdam.Subroutines.output import generate_output_files
from tests.make_atom_rec import make_exp_mmcif_rec

class TestClass(unittest.TestCase):

    def test_download_mmcif(self):
        """
        Checks that downloads mmCIF files correctly, and raises the appropriate
        errors e.g. if the file does not exist
        """

        import os
        import shutil

        exp_results = {'2BN1': False,
                       '1ABC': True}
        for pdb, exp_exit in exp_results.items():
            if os.path.isdir('tests/temp_files/'):
                shutil.rmtree('tests/temp_files/')
            act_exit = download_mmcif(
                pdb, 'tests/temp_files/', 'tests/temp_files/test.cif'
            )
            self.assertEqual(exp_exit, act_exit)

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

        # Expect CRYST1 line to be constructed correctly
        exp_cryst1_line_1 = 'CRYST1   77.900   77.900   77.900  90.00  90.00  90.00 I 21 3       24          \n'
        act_cryst1_line_1, exit_1 = make_cryst1_line_from_mmcif(space_group_info_1, False)
        self.assertEqual(exp_cryst1_line_1, act_cryst1_line_1)
        self.assertFalse(exit_1)  

        # Expect error because dimension is too large (should be 9 characters,
        # with 3 characters reserved for after the decimal point)
        space_group_info_2 = copy.deepcopy(space_group_info_1)
        space_group_info_2[1] = '_cell.length_a       123477.900'
        exp_cryst1_line_2 = ''
        act_cryst1_line_2, exit_2 = make_cryst1_line_from_mmcif(space_group_info_2, False)
        self.assertEqual(exp_cryst1_line_2, act_cryst1_line_2)
        self.assertTrue(exit_2)

        # Expect error because dimension is not a numeric value
        space_group_info_3 = copy.deepcopy(space_group_info_1)
        space_group_info_3[2] = '_cell.length_b                A'
        exp_cryst1_line_3 = ''
        act_cryst1_line_3, exit_3 = make_cryst1_line_from_mmcif(space_group_info_3, False)
        self.assertEqual(exp_cryst1_line_3, act_cryst1_line_3)
        self.assertTrue(exit_3)

        # Expect CRYST1 line to be constructed correctly
        space_group_info_4 = copy.deepcopy(space_group_info_1)
        space_group_info_4[3] = '_cell.length_c 18.5'
        exp_cryst1_line_4 = 'CRYST1   77.900   77.900   18.500  90.00  90.00  90.00 I 21 3       24          \n'
        act_cryst1_line_4, exit_4 = make_cryst1_line_from_mmcif(space_group_info_4, False)
        self.assertEqual(exp_cryst1_line_4, act_cryst1_line_4)
        self.assertFalse(exit_4)

        # Expect CRYST1 line to be constructed correctly
        space_group_info_5 = copy.deepcopy(space_group_info_1)
        space_group_info_5[4] = '_cell.angle_alpha        15.00'
        exp_cryst1_line_5 = 'CRYST1   77.900   77.900   77.900  15.00  90.00  90.00 I 21 3       24          \n'
        act_cryst1_line_5, exit_5 = make_cryst1_line_from_mmcif(space_group_info_5, False)
        self.assertEqual(exp_cryst1_line_5, act_cryst1_line_5)
        self.assertFalse(exit_5)

        # Expect CRYST1 line to be constructed correctly
        space_group_info_6 = copy.deepcopy(space_group_info_1)
        space_group_info_6[5] = '_cell.angle_beta         15.5678'
        exp_cryst1_line_6 = 'CRYST1   77.900   77.900   77.900  90.00  15.57  90.00 I 21 3       24          \n'
        act_cryst1_line_6, exit_6 = make_cryst1_line_from_mmcif(space_group_info_6, False)
        self.assertEqual(exp_cryst1_line_6, act_cryst1_line_6)
        self.assertFalse(exit_6)

        # Expect error because dimension is too large (should be 7 characters,
        # with 2 characters reserved for after the decimal point)
        space_group_info_7 = copy.deepcopy(space_group_info_1)
        space_group_info_7[5] = '_cell.angle_gamma         123456'
        exp_cryst1_line_7 = ''
        act_cryst1_line_7, exit_7 = make_cryst1_line_from_mmcif(space_group_info_7, False)
        self.assertEqual(exp_cryst1_line_7, act_cryst1_line_7)
        self.assertTrue(exit_7)

        # Expect error because H-M space group name is longer than 11 characters
        space_group_info_8 = copy.deepcopy(space_group_info_1)
        space_group_info_8[7] = '_symmetry.space_group_name_H-M \'123456789101\''
        exp_cryst1_line_8 = ''
        act_cryst1_line_8, exit_8 = make_cryst1_line_from_mmcif(space_group_info_8, False)
        self.assertEqual(exp_cryst1_line_8, act_cryst1_line_8)
        self.assertTrue(exit_8)

        # Expect error because Z_PDB (= the number of polymeric chains in a
        # unit cell) is longer than 4 characters
        space_group_info_9 = copy.deepcopy(space_group_info_1)
        space_group_info_9[8] = '_cell.Z_PDB 12345'
        exp_cryst1_line_9 = ''
        act_cryst1_line_9, exit_9 = make_cryst1_line_from_mmcif(space_group_info_9, False)
        self.assertEqual(exp_cryst1_line_9, act_cryst1_line_9)
        self.assertTrue(exit_9)

        # Expect error because one of the expected unit cell dimensions is not
        # supplied
        space_group_info_10 = copy.deepcopy(space_group_info_1)
        space_group_info_10[1] = ''
        exp_cryst1_line_10 = ''
        act_cryst1_line_10, exit_10 = make_cryst1_line_from_mmcif(space_group_info_10, False)
        self.assertEqual(exp_cryst1_line_10, act_cryst1_line_10)
        self.assertTrue(exit_10)

    def test_parse_atomrec_from_mmcif(self):
        """
        Checks that ATOM/HETATM records in mmCIF format files are parsed as
        expected
        """

        import os
        import shutil

        exp_atom_rec_obj, exp_atom_rec_str = make_exp_mmcif_rec()
        if not os.path.isdir('tests/temp_files/'):
            os.mkdir('tests/temp_files/')
        with open('tests/temp_files/input.cif', 'w') as f:
            f.write('data_1ABC\n#\nloop_\n')
            f.write('\n'.join(exp_atom_rec_str))
            f.write('\n#')

        act_atom_rec, exit = parse_atom_rec_from_mmcif(
            exp_atom_rec_str, 'tests/temp_files/input.cif', False
        )
        for index in range(len(act_atom_rec)):
            self.assertEqual(exp_atom_rec_obj[index], act_atom_rec[index])
        self.assertFalse(exit)

        shutil.rmtree('tests/temp_files/')

    def test_write_mmcif_file(self):
        """
        Tests that atom objects are correctly converted into mmCIF format
        """

        import os
        import pandas as pd
        import shutil

        if not os.path.isdir('tests/temp_files/'):
            os.mkdir('tests/temp_files/')

        atoms_list = []
        # Taken from 2BN1 and 2QTZ
        new_atom_1 = atom()
        new_atom_1.lineID = 'ATOM'
        new_atom_1.atomNum = 2
        new_atom_1.atomType = 'CA'
        new_atom_1.conformer = ''
        new_atom_1.resiType = 'GLY'
        new_atom_1.chainID = 'A'
        new_atom_1.resiNum = 1
        new_atom_1.insCode = ''
        new_atom_1.xyzCoords = [[13.602], [45.768], [30.728]]
        new_atom_1.occupancy = 1.0
        new_atom_1.bFactor = 19.08
        new_atom_1.element = 'C'
        new_atom_1.charge = ''
        new_atom_1.origResiNum = '1'
        new_atom_1.origResiType = 'GLY'
        new_atom_1.origChainID = 'A'
        new_atom_1.origAtomType = 'CA'
        new_atom_1.pd = 19
        new_atom_1.avrg_bf = 26.2
        new_atom_1.bd = 1.3
        new_atom_1.protein = False
        new_atom_1.na = False
        new_atom_1.chain_len = 64
        atoms_list.append(new_atom_1)

        new_atom_2 = atom()
        new_atom_2.lineID = 'HETATM'
        new_atom_2.atomNum = 438
        new_atom_2.atomType = 'O'
        new_atom_2.conformer = ''
        new_atom_2.resiType = 'HOH'
        new_atom_2.chainID = 'B'
        new_atom_2.resiNum = 2001
        new_atom_2.insCode = ''
        new_atom_2.xyzCoords = [[13.128], [46.298], [34.342]]
        new_atom_2.occupancy = 1.0
        new_atom_2.bFactor = 33.1
        new_atom_2.element = 'O'
        new_atom_2.charge = ''
        new_atom_2.origResiNum = ''
        new_atom_2.origResiType = 'HOH'
        new_atom_2.origChainID = 'C'
        new_atom_2.origAtomType = 'O'
        new_atom_2.pd = 30
        new_atom_2.avrg_bf = 25.6
        new_atom_2.bd = 0.97
        new_atom_2.protein = True
        new_atom_2.na = True
        new_atom_2.chain_len = 0
        atoms_list.append(new_atom_2)

        new_atom_3 = atom()
        new_atom_3.lineID = 'HETATM'
        new_atom_3.atomNum = 3546
        new_atom_3.atomType = 'PA'
        new_atom_3.conformer = 'X'
        new_atom_3.resiType = 'FAD'
        new_atom_3.chainID = 'B'
        new_atom_3.resiNum = 700
        new_atom_3.insCode = ''
        new_atom_3.xyzCoords = [[29.643], [29.7], [51.402]]
        new_atom_3.occupancy = 1.0
        new_atom_3.bFactor = 20.52
        new_atom_3.element = 'P'
        new_atom_3.charge = ''
        new_atom_3.origResiNum = ''
        new_atom_3.origResiType = 'FAD'
        new_atom_3.origChainID = 'B'
        new_atom_3.origAtomType = 'PA'
        new_atom_3.pd = 26
        new_atom_3.avrg_bf = 20.5
        new_atom_3.bd = 2.0
        new_atom_3.protein = False
        new_atom_3.na = True
        new_atom_3.chain_len = 1
        atoms_list.append(new_atom_3)

        exp_cif_lines = [
            '#',
            'loop_',
            '_atom_site.group_PDB',
            '_atom_site.id',
            '_atom_site.type_symbol',
            '_atom_site.label_atom_id',
            '_atom_site.label_alt_id',
            '_atom_site.label_comp_id',
            '_atom_site.label_asym_id',
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
            'ATOM   2    C CA . GLY A 1 ? 13.602 45.768 30.728 1.0 19.08 1.3  19 26.2 ? 1    GLY A CA',
            'HETATM 438  O O  . HOH C . ? 13.128 46.298 34.342 1.0 33.1  0.97 30 25.6 ? 2001 HOH B O ',
            'HETATM 3546 P PA X FAD B . ? 29.643 29.7   51.402 1.0 20.52 2.0  26 20.5 ? 700  FAD B PA',
            '#',
            ''
        ]

        empty_df = pd.DataFrame({})
        output = generate_output_files('tests/temp_files/test', '', empty_df)
        output.write_output_cif(atoms_list)
        with open('tests/temp_files/test_BDamage.cif', 'r') as f:
            act_cif_lines = f.read().split('\n')

        self.assertEqual(exp_cif_lines, act_cif_lines)

        shutil.rmtree('tests/temp_files/')
