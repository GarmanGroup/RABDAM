
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

# An outer layer to the pipeline scripts. Depending upon the flags specified
# in the command line input, this script will run either the complete / a
# subsection of the pipeline.

# python -m unittest tests/test_atom_filtering.py

import unittest
from rabdam.Subroutines.parsePDB import atom, b_damage_atom_list
from rabdam.Subroutines.PDBCUR import clean_atom_rec

class TestClass(unittest.TestCase):

    def make_atoms_list(self):
        """
        """

        atoms_list = []
        # Taken from 2BN1 and 2QTZ
        new_atom_1 = atom()
        new_atom_1.lineID = 'ATOM'
        new_atom_1.atomNum = 1
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
        new_atom_1.origResiNum = ''
        new_atom_1.origResiType = 'GLY'
        new_atom_1.origChainID = 'A'
        new_atom_1.origAtomType = 'CA'
        new_atom_1.pd = 19
        new_atom_1.avrg_bf = 26.2
        new_atom_1.bd = 1.3
        atoms_list.append(new_atom_1)

        new_atom_2 = atom()
        new_atom_2.lineID = 'ATOM'
        new_atom_2.atomNum = 2
        new_atom_2.atomType = 'CA'
        new_atom_2.conformer = ''
        new_atom_2.resiType = 'GLY'
        new_atom_2.chainID = 'A'
        new_atom_2.resiNum = 2
        new_atom_2.insCode = ''
        new_atom_2.xyzCoords = [[13.128], [46.298], [34.342]]
        new_atom_2.occupancy = 1.0
        new_atom_2.bFactor = 33.1
        new_atom_2.element = 'C'
        new_atom_2.charge = ''
        new_atom_2.origResiNum = ''
        new_atom_2.origResiType = 'GLY'
        new_atom_2.origChainID = 'A'
        new_atom_2.origAtomType = 'CA'
        new_atom_2.pd = 30
        new_atom_2.avrg_bf = 25.6
        new_atom_2.bd = 0.97
        atoms_list.append(new_atom_2)

        new_atom_3 = atom()
        new_atom_3.lineID = 'HETATM'
        new_atom_3.atomNum = 3
        new_atom_3.atomType = 'PA'
        new_atom_3.conformer = 'X'
        new_atom_3.resiType = 'FAD'
        new_atom_3.chainID = 'B'
        new_atom_3.resiNum = 2
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
        atoms_list.append(new_atom_3)

        return atoms_list

    def test_clean_atom_records(self):
        """
        Checks that ATOM/HETATM records are filtered appropriately to remove
        hydrogen atoms, 0 occupancy atoms, retain only a single conformer
        per-residue, and check that disulfide bonds have not been refined to
        sub-1 occupancy
        """

        import copy
        import os
        import shutil

        if not os.path.isdir('tests/temp_files/'):
            os.mkdir('tests/temp_files/')

        atoms_list_1 = self.make_atoms_list()
        exp_atoms_list_1 = copy.deepcopy(atoms_list_1)
        exp_file_content_1 = [
            '',
            'ATOM      1  CA  GLY A   1      13.602  45.768  30.728  1.00 19.08           C  ',
            'ATOM      2  CA  GLY A   2      13.128  46.298  34.342  1.00 33.10           C  ',
            'TER                                                                             ',
            'HETATM    3  PA XFAD B   2      29.643  29.700  51.402  1.00 20.52           P  ',
            'TER                                                                             ',
            ''
        ]
        exit_1, pause_1, act_atoms_list_1, file_1 = clean_atom_rec(
            atoms_list_1, {}, ['GLY', 'FAD'], '', 'tests/temp_files/test'
        )
        self.assertFalse(exit_1)
        self.assertFalse(pause_1)
        self.assertEqual(exp_atoms_list_1, act_atoms_list_1)
        with open('tests/temp_files/test_asymmetric_unit.pdb', 'r') as f:
            act_file_content_1 = f.read().split('\n')
            self.assertEqual(exp_file_content_1, act_file_content_1)

        # Hydrogen atoms
        atoms_list_2 = copy.deepcopy(atoms_list_1)
        atoms_list_2[0].element = 'H'
        exp_atoms_list_2 = copy.deepcopy(atoms_list_2)[1:]
        exp_file_content_2 = [
            '',
            'ATOM      2  CA  GLY A   2      13.128  46.298  34.342  1.00 33.10           C  ',
            'TER                                                                             ',
            'HETATM    3  PA XFAD B   2      29.643  29.700  51.402  1.00 20.52           P  ',
            'TER                                                                             ',
            ''
        ]
        exit_2, pause_2, act_atoms_list_2, file_2 = clean_atom_rec(
            atoms_list_2, {}, ['GLY', 'FAD'], '', 'tests/temp_files/test'
        )
        self.assertFalse(exit_2)
        self.assertFalse(pause_2)
        self.assertEqual(exp_atoms_list_2, act_atoms_list_2)
        with open('tests/temp_files/test_asymmetric_unit.pdb', 'r') as f:
            act_file_content_2 = f.read().split('\n')
            self.assertEqual(exp_file_content_2, act_file_content_2)

        # 0 occupancy
        atoms_list_3 = copy.deepcopy(atoms_list_1)
        atoms_list_3[1].occupancy = 0
        exp_atoms_list_3 = [copy.deepcopy(atoms_list_3)[0], copy.deepcopy(atoms_list_3)[2]]
        exp_file_content_3 = [
            '',
            'ATOM      1  CA  GLY A   1      13.602  45.768  30.728  1.00 19.08           C  ',
            'TER                                                                             ',
            'HETATM    3  PA XFAD B   2      29.643  29.700  51.402  1.00 20.52           P  ',
            'TER                                                                             ',
            ''
        ]
        exit_3, pause_3, act_atoms_list_3, file_3 = clean_atom_rec(
            atoms_list_3, {}, ['GLY', 'FAD'], '', 'tests/temp_files/test'
        )
        self.assertFalse(exit_3)
        self.assertFalse(pause_3)
        self.assertEqual(exp_atoms_list_3, act_atoms_list_3)
        with open('tests/temp_files/test_asymmetric_unit.pdb', 'r') as f:
            act_file_content_3 = f.read().split('\n')
            self.assertEqual(exp_file_content_3, act_file_content_3)

        # Negative Bfactor
        atoms_list_4 = copy.deepcopy(atoms_list_1)
        atoms_list_4[2].bFactor = -0.1
        exp_atoms_list_4 = copy.deepcopy(atoms_list_4)[:2]
        exp_file_content_4 = [
            '',
            'ATOM      1  CA  GLY A   1      13.602  45.768  30.728  1.00 19.08           C  ',
            'ATOM      2  CA  GLY A   2      13.128  46.298  34.342  1.00 33.10           C  ',
            'TER                                                                             ',
            ''
        ]
        exit_4, pause_4, act_atoms_list_4, file_4 = clean_atom_rec(
            atoms_list_4, {}, ['GLY', 'FAD'], '', 'tests/temp_files/test'
        )
        self.assertFalse(exit_4)
        self.assertFalse(pause_4)
        self.assertEqual(exp_atoms_list_4, act_atoms_list_4)
        with open('tests/temp_files/test_asymmetric_unit.pdb', 'r') as f:
            act_file_content_4 = f.read().split('\n')
            self.assertEqual(exp_file_content_4, act_file_content_4)

        # Single sub-1 occupancy
        atoms_list_5 = copy.deepcopy(atoms_list_1)
        atoms_list_5[2].occupancy = 0.7
        exp_atoms_list_5 = copy.deepcopy(atoms_list_5)
        exp_file_content_5 = [
            '',
            'ATOM      1  CA  GLY A   1      13.602  45.768  30.728  1.00 19.08           C  ',
            'ATOM      2  CA  GLY A   2      13.128  46.298  34.342  1.00 33.10           C  ',
            'TER                                                                             ',
            'HETATM    3  PA XFAD B   2      29.643  29.700  51.402  0.70 20.52           P  ',
            'TER                                                                             ',
            ''
        ]
        exit_5, pause_5, act_atoms_list_5, file_5 = clean_atom_rec(
            atoms_list_5, {}, ['GLY', 'FAD'], '', 'tests/temp_files/test'
        )
        self.assertFalse(exit_5)
        self.assertTrue(pause_5)
        self.assertEqual(exp_atoms_list_5, act_atoms_list_5)
        with open('tests/temp_files/test_asymmetric_unit.pdb', 'r') as f:
            act_file_content_5 = f.read().split('\n')
            self.assertEqual(exp_file_content_5, act_file_content_5)

        # Alternate conformers that sum to 1
        atoms_list_6 = copy.deepcopy(atoms_list_1)
        atoms_list_6[0].conformer = 'A'
        atoms_list_6[0].occupancy = 0.25
        atoms_list_6[1].resiNum = 1
        atoms_list_6[1].conformer = 'B'
        atoms_list_6[1].occupancy = 0.75
        exp_atoms_list_6 = copy.deepcopy(atoms_list_6)[1:]
        exp_file_content_6 = [
            '',
            'ATOM      2  CA BGLY A   1      13.128  46.298  34.342  0.75 33.10           C  ',
            'TER                                                                             ',
            'HETATM    3  PA XFAD B   2      29.643  29.700  51.402  1.00 20.52           P  ',
            'TER                                                                             ',
            ''
        ]
        exit_6, pause_6, act_atoms_list_6, file_6 = clean_atom_rec(
            atoms_list_6, {}, ['GLY', 'FAD'], '', 'tests/temp_files/test'
        )
        self.assertFalse(exit_6)
        self.assertFalse(pause_6)
        self.assertEqual(exp_atoms_list_6, act_atoms_list_6)
        with open('tests/temp_files/test_asymmetric_unit.pdb', 'r') as f:
            act_file_content_6 = f.read().split('\n')
            self.assertEqual(exp_file_content_6, act_file_content_6)

        # Alternate conformers that don't sum to 1
        atoms_list_7 = copy.deepcopy(atoms_list_1)
        atoms_list_7[0].conformer = 'A'
        atoms_list_7[0].occupancy = 0.25
        atoms_list_7[1].resiNum = 1
        atoms_list_7[1].conformer = 'B'
        atoms_list_7[1].occupancy = 0.74
        exp_atoms_list_7 = copy.deepcopy(atoms_list_7)[1:]
        exp_file_content_7 = [
            '',
            'ATOM      2  CA BGLY A   1      13.128  46.298  34.342  0.74 33.10           C  ',
            'TER                                                                             ',
            'HETATM    3  PA XFAD B   2      29.643  29.700  51.402  1.00 20.52           P  ',
            'TER                                                                             ',
            ''
        ]
        exit_7, pause_7, act_atoms_list_7, file_7 = clean_atom_rec(
            atoms_list_7, {}, ['GLY', 'FAD'], '', 'tests/temp_files/test'
        )
        self.assertFalse(exit_7)
        self.assertTrue(pause_7)
        self.assertEqual(exp_atoms_list_7, act_atoms_list_7)
        with open('tests/temp_files/test_asymmetric_unit.pdb', 'r') as f:
            act_file_content_7 = f.read().split('\n')
            self.assertEqual(exp_file_content_7, act_file_content_7)

        # Disulfide bonds
        atoms_list_8 = copy.deepcopy(atoms_list_1)
        atoms_list_8[0].resiType = 'CYS'
        atoms_list_8[0].occupancy = 0.25
        exp_atoms_list_8 = copy.deepcopy(atoms_list_8)
        exp_file_content_8 = [
            '',
            'ATOM      1  CA  CYS A   1      13.602  45.768  30.728  0.25 19.08           C  ',
            'ATOM      2  CA  GLY A   2      13.128  46.298  34.342  1.00 33.10           C  ',
            'TER                                                                             ',
            'HETATM    3  PA XFAD B   2      29.643  29.700  51.402  1.00 20.52           P  ',
            'TER                                                                             ',
            ''
        ]
        exit_8, pause_8, act_atoms_list_8, file_8 = clean_atom_rec(
            atoms_list_8, {1: [['A', 1, ''], ['B', 7, '']]}, ['GLY', 'FAD'],
            '', 'tests/temp_files/test'
        )
        self.assertFalse(exit_8)
        self.assertTrue(pause_8)
        self.assertEqual(exp_atoms_list_8, act_atoms_list_8)
        with open('tests/temp_files/test_asymmetric_unit.pdb', 'r') as f:
            act_file_content_8 = f.read().split('\n')
            self.assertEqual(exp_file_content_8, act_file_content_8)

        # No atoms left after filtering
        atoms_list_9 = copy.deepcopy(atoms_list_1)
        atoms_list_9[0].element = 'H'
        atoms_list_9[1].occupancy = 1.4
        atoms_list_9[2].bFactor = 0
        exp_atoms_list_9 = []
        exp_file_content_9 = ['', '']
        exit_9, pause_9, act_atoms_list_9, file_9 = clean_atom_rec(
            atoms_list_9, {}, ['GLY', 'FAD'], '', 'tests/temp_files/test'
        )
        self.assertTrue(exit_9)
        self.assertFalse(pause_9)
        self.assertEqual(exp_atoms_list_9, act_atoms_list_9)
        with open('tests/temp_files/test_asymmetric_unit.pdb', 'r') as f:
            act_file_content_9 = f.read().split('\n')
            self.assertEqual(exp_file_content_9, act_file_content_9)

        shutil.rmtree('tests/temp_files/')

    def test_create_list_of_atoms_for_bdamage_analysis(self):
        """
        Checks that input list of atom objects is correctly filtered
        according to user-specified program inputs
        """

        import copy

        atoms_list_1 = self.make_atoms_list()

        # Remove/keep HETATM
        atoms_list_1[2].resiType = 'MSE'
        act_atoms_list_1 = b_damage_atom_list(
            atoms_list_1, [], False, 'protein', [], []
        )
        exp_atoms_list_1 = copy.deepcopy(atoms_list_1)[:2]
        self.assertEqual(exp_atoms_list_1, act_atoms_list_1)

        atoms_list_2 = copy.deepcopy(atoms_list_1)
        atoms_list_2[2].resiType = 'MSE'
        act_atoms_list_2 = b_damage_atom_list(
            atoms_list_2, ['GLY', 'MSE'], False, 'protein', [], []
        )
        exp_atoms_list_2 = copy.deepcopy(atoms_list_2)
        self.assertEqual(exp_atoms_list_2, act_atoms_list_2)

        atoms_list_3 = copy.deepcopy(atoms_list_1)
        atoms_list_3[0].lineID = 'HETATM'
        act_atoms_list_3 = b_damage_atom_list(
            atoms_list_3, ['GLY', 'MSE'], False, 'protein', [], []
        )
        exp_atoms_list_3 = copy.deepcopy(atoms_list_3)[1:]
        self.assertEqual(exp_atoms_list_3, act_atoms_list_3)

        # Keep protein/NA
        atoms_list_4 = copy.deepcopy(atoms_list_1)
        atoms_list_4[1].resiType = 'DU'
        atoms_list_4[2].resiType = '8OG'
        act_atoms_list_4 = b_damage_atom_list(
            atoms_list_4, ['GLY', 'DU', '8OG'], False, 'protein', [], []
        )
        exp_atoms_list_4 = copy.deepcopy(atoms_list_4)[0:1]
        self.assertEqual(exp_atoms_list_4, act_atoms_list_4)

        atoms_list_5 = copy.deepcopy(atoms_list_1)
        atoms_list_5[1].resiType = 'DU'
        atoms_list_5[2].resiType = '8OG'
        act_atoms_list_5 = b_damage_atom_list(
            atoms_list_5, ['GLY', 'DU'], False, 'na', [], []
        )
        exp_atoms_list_5 = copy.deepcopy(atoms_list_5)[1:2]
        self.assertEqual(exp_atoms_list_5, act_atoms_list_5)

        atoms_list_6 = copy.deepcopy(atoms_list_1)
        atoms_list_6[1].resiType = 'DU'
        atoms_list_6[2].resiType = '8OG'
        act_atoms_list_6 = b_damage_atom_list(
            atoms_list_6, ['GLY', 'DU'], True, 'na', [], []
        )
        exp_atoms_list_6 = copy.deepcopy(atoms_list_6)[1:]
        self.assertEqual(exp_atoms_list_6, act_atoms_list_6)

        # Remove atoms
        atoms_list_7 = copy.deepcopy(atoms_list_1)
        act_atoms_list_7 = b_damage_atom_list(
            atoms_list_7, [], True, 'na', [], ['1', '3']
        )
        exp_atoms_list_7 = copy.deepcopy(atoms_list_7)[1:2]
        self.assertEqual(exp_atoms_list_7, act_atoms_list_7)

        atoms_list_8 = copy.deepcopy(atoms_list_1)
        act_atoms_list_8 = b_damage_atom_list(
            atoms_list_8, [], True, 'na', [], ['GLY']
        )
        exp_atoms_list_8 = copy.deepcopy(atoms_list_8)[2:]
        self.assertEqual(exp_atoms_list_8, act_atoms_list_8)

        # Add atoms
        atoms_list_9 = copy.deepcopy(atoms_list_1)
        act_atoms_list_9 = b_damage_atom_list(
            atoms_list_9, [], True, 'na', ['2'], ['GLY']
        )
        exp_atoms_list_9 = copy.deepcopy(atoms_list_9)[1:]
        self.assertEqual(exp_atoms_list_9, act_atoms_list_9)
