
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

# python -m unittest tests/test_pdb_parsing.py

import unittest
from rabdam.Subroutines.parsePDB import atom
from rabdam.Subroutines.PDBCUR import (
    find_disulfides_from_pdb, parse_seqres_from_pdb, parse_atom_rec_from_pdb
)
from rabdam.Subroutines.makeDataFrame import makePDB

class TestClass(unittest.TestCase):

    def test_find_disulfides_from_pdb(self):
        """
        Checks that disulfide bond records in PDB format files are parsed as
        expected
        """

        import copy

        # Taken from 2BN1
        disulfide_rec_1 = [
            'SSBOND   1 CYS A    6    CYS A   11                          1555   1555  2.16  ',
            'SSBOND   2 CYS A    7    CYS B    7                          1555   1555  2.13  ',
            'SSBOND   3 CYS A   20    CYS B   19                          1555   1555  2.02  '
        ]
        exp_disulfides_1 = {1: [['A', 6, ''], ['A', 11, '']],
                            2: [['A', 7, ''], ['B', 7, '']],
                            3: [['A', 20, ''], ['B', 19, '']]}
        act_disulfides_1, exit_1 = find_disulfides_from_pdb(disulfide_rec_1, False)
        self.assertDictEqual(act_disulfides_1, exp_disulfides_1)
        self.assertFalse(exit_1)

        act_disulfides_2, exit_2 = find_disulfides_from_pdb(disulfide_rec_1, True)
        self.assertDictEqual(act_disulfides_2, exp_disulfides_1)
        self.assertTrue(exit_2)


        disulfide_rec_3 = copy.deepcopy(disulfide_rec_1)
        disulfide_rec_3[2] = 'SSBOND   3 CYS A    X    CYS B   19                          1555   1555  2.02  '
        act_disulfides_3, exit_3 = find_disulfides_from_pdb(disulfide_rec_3, False)
        self.assertDictEqual(act_disulfides_3, {})
        self.assertTrue(exit_3)

        disulfide_rec_4 = copy.deepcopy(disulfide_rec_1)
        disulfide_rec_4[2] = 'SSBOND'
        act_disulfides_4, exit_4 = find_disulfides_from_pdb(disulfide_rec_4, False)
        self.assertDictEqual(act_disulfides_4, {})
        self.assertTrue(exit_4)

    def test_parse_seqres_from_pdb(self):
        """
        Checks that SEQRES records in PDB format files are parsed as expected
        """

        import copy

        # Taken from 2BN1
        seqres_rec_1 = [
            'SEQRES   1 A   21  GLY ILE VAL GLU GLN CYS CYS ALA SER VAL CYS SER LEU          ',
            'SEQRES   2 A   21  TYR GLN LEU GLU ASN TYR CYS ASN                              ',
            'SEQRES   1 B   30  PHE VAL ASN GLN HIS LEU CYS GLY SER HIS LEU VAL GLU          ',
            'SEQRES   2 B   30  ALA LEU TYR LEU VAL CYS GLY GLU ARG GLY PHE PHE TYR          ',
            'SEQRES   3 B   30  THR PRO LYS ALA                                              '
        ]
        exp_seq_1 = [
            'GLY', 'ILE', 'VAL', 'GLU', 'GLN', 'CYS', 'ALA', 'SER', 'LEU',
            'TYR', 'ASN', 'PHE', 'HIS', 'ARG', 'THR', 'PRO', 'LYS'
        ]
        act_seq_1, exit_1 = parse_seqres_from_pdb(seqres_rec_1, False)
        self.assertEqual(exp_seq_1, act_seq_1)
        self.assertFalse(exit_1)

        act_seq_2, exit_2 = parse_seqres_from_pdb(seqres_rec_1, True)
        self.assertEqual(exp_seq_1, act_seq_2)
        self.assertTrue(exit_2)

        seqres_rec_3 = copy.deepcopy(seqres_rec_1)
        seqres_rec_3[1] = 'SEQRES   2 A   21  TYR GLN LEU GLU ASN TYR CYS ASN'
        act_seq_3, exit_3 = parse_seqres_from_pdb(seqres_rec_3, False)
        self.assertEqual(exp_seq_1, act_seq_3)
        self.assertFalse(exit_3)

        seqres_rec_4 = copy.deepcopy(seqres_rec_1)
        seqres_rec_4[1] = 'SEQRES   2 A   21'
        exp_seq_4 = []
        act_seq_4, exit_4 = parse_seqres_from_pdb(seqres_rec_4, False)
        self.assertEqual(exp_seq_4, act_seq_4)
        self.assertTrue(exit_4)

    def test_parse_atomrec_from_pdb(self):
        """
        Checks that PDB format files are parsed as expected
        """

        import numpy as np

        # Taken from 2BN1
        atom_lines_1 = ['ATOM      2  CA  GLY A   1      13.602  45.768  30.728  1.00 19.08           C  ']

        exp_atom_rec_1 = []
        exp_atom_1 = atom()
        exp_atom_1.lineID = 'ATOM'
        exp_atom_1.atomNum = 2
        exp_atom_1.atomType = 'CA'
        exp_atom_1.conformer = ''
        exp_atom_1.resiType = 'GLY'
        exp_atom_1.chainID = 'A'
        exp_atom_1.resiNum = 1
        exp_atom_1.insCode = ''
        exp_atom_1.xyzCoords = [[13.602], [45.768], [30.728]]
        exp_atom_1.occupancy = 1.0
        exp_atom_1.bFactor = 19.08
        exp_atom_1.element = 'C'
        exp_atom_1.charge = ''
        exp_atom_1.origResiNum = '1'
        exp_atom_1.origResiType = 'GLY'
        exp_atom_1.origChainID = 'A'
        exp_atom_1.origAtomType = 'CA'
        exp_atom_1.pd = None
        exp_atom_1.avrg_bf = None
        exp_atom_1.bd = None
        exp_atom_rec_1.append(exp_atom_1)

        act_atom_rec_1, exit_1 = parse_atom_rec_from_pdb(atom_lines_1, False)
        for index in range(len(act_atom_rec_1)):
            self.assertEqual(exp_atom_rec_1[index], act_atom_rec_1[index])
        self.assertFalse(exit_1)

        act_atom_rec_2, exit_2 = parse_atom_rec_from_pdb(atom_lines_1, True)
        for index in range(len(act_atom_rec_2)):
            self.assertEqual(exp_atom_rec_1[index], act_atom_rec_2[index])
        self.assertTrue(exit_2)

        atom_lines_3 = ['ATOM      X  CA  GLY A   1      13.602  45.768  30.728  1.00 19.08           C  ']
        exp_atom_rec_3 = []
        act_atom_rec_3, exit_3 = parse_atom_rec_from_pdb(atom_lines_3, False)
        self.assertEqual(exp_atom_rec_3, act_atom_rec_3)
        self.assertTrue(exit_3)

        atom_lines_4 = ['ATOM      2  CA  GLY A   X      13.602  45.768  30.728  1.00 19.08           C  ']
        exp_atom_rec_4 = []
        act_atom_rec_4, exit_4 = parse_atom_rec_from_pdb(atom_lines_4, False)
        self.assertEqual(exp_atom_rec_4, act_atom_rec_4)
        self.assertTrue(exit_4)

        atom_lines_5 = ['ATOM      2  CA  GLY A   1           X  45.768  30.728  1.00 19.08           C  ']
        exp_atom_rec_5 = []
        act_atom_rec_5, exit_5 = parse_atom_rec_from_pdb(atom_lines_5, False)
        self.assertEqual(exp_atom_rec_5, act_atom_rec_5)
        self.assertTrue(exit_5)

        atom_lines_6 = ['ATOM      2  CA  GLY A   1      13.602  45.768  30.728     X 19.08           C  ']
        exp_atom_rec_6 = []
        act_atom_rec_6, exit_6 = parse_atom_rec_from_pdb(atom_lines_6, False)
        self.assertEqual(exp_atom_rec_6, act_atom_rec_6)
        self.assertTrue(exit_6)

        atom_lines_7 = ['ATOM      2  CA  GLY A   1      13.602  45.768  30.728  1.00     X           C  ']
        exp_atom_rec_7 = []
        act_atom_rec_7, exit_7 = parse_atom_rec_from_pdb(atom_lines_7, False)
        self.assertEqual(exp_atom_rec_7, act_atom_rec_7)
        self.assertTrue(exit_7)

    def test_write_pdb_file(self):
        """
        Tests that atom objects are correctly converted into PDB format, and TER
        cards are inserted at appropriate locations
        """

        import os
        import copy
        import shutil

        if not os.path.isdir('tests/temp_files/'):
            os.mkdir('tests/temp_files/')

        atoms_list_1 = []

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
        new_atom_1.pd = None
        new_atom_1.avrg_bf = None
        new_atom_1.bd = 19.08
        atoms_list_1.append(new_atom_1)

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
        new_atom_2.pd = None
        new_atom_2.avrg_bf = None
        new_atom_2.bd = 33.1
        atoms_list_1.append(new_atom_2)

        new_atom_3 = atom()
        new_atom_3.lineID = 'HETATM'
        new_atom_3.atomNum = 3546
        new_atom_3.atomType = 'PA'
        new_atom_3.conformer = ''
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
        new_atom_3.pd = None
        new_atom_3.avrg_bf = None
        new_atom_3.bd = 20.52
        atoms_list_1.append(new_atom_3)

        exit_1 = makePDB(
            ['Header'], atoms_list_1, ['Footer'], [], 'tests/temp_files/test.pdb', 'bfactor'
        )
        with open('tests/temp_files/test.pdb', 'r') as f:
            act_pdb_lines_1 = f.read().split('\n')
        exp_pdb_lines_1 = [
            'Header',
            'ATOM      2  CA  GLY A   1      13.602  45.768  30.728  1.00 19.08           C  ',
            'HETATM  438  O   HOH B2001      13.128  46.298  34.342  1.00 33.10           O  ',
            'HETATM 3546  PA  FAD B 700      29.643  29.700  51.402  1.00 20.52           P  ',
            'Footer',
            ''
        ]
        self.assertEqual(act_pdb_lines_1, exp_pdb_lines_1)
        self.assertFalse(exit_1)

        atoms_list_2 = copy.deepcopy(atoms_list_1)
        atoms_list_2[0].lineID = 'HETATMM'
        exit_2 = makePDB(
            ['Header'], atoms_list_2, ['Footer'], [], 'tests/temp_files/test.pdb', 'bfactor'
        )
        with open('tests/temp_files/test.pdb', 'r') as f:
            act_pdb_lines_2 = f.read().split('\n')
        exp_pdb_lines_2 = ['Header', 'Footer', '']
        self.assertEqual(act_pdb_lines_2, exp_pdb_lines_2)
        self.assertTrue(exit_2)

        atoms_list_3 = copy.deepcopy(atoms_list_1)
        atoms_list_3[0].atomType = 'CG11'
        exit_3 = makePDB(
            ['Header'], atoms_list_3, ['Footer'], [], 'tests/temp_files/test.pdb', 'bfactor'
        )
        with open('tests/temp_files/test.pdb', 'r') as f:
            act_pdb_lines_3 = f.read().split('\n')
        exp_pdb_lines_3 = ['Header', 'Footer', '']
        self.assertEqual(act_pdb_lines_3, exp_pdb_lines_3)
        self.assertTrue(exit_3)

        atoms_list_4 = copy.deepcopy(atoms_list_1)
        atoms_list_4[0].conformer = 'XYZ'
        exit_4 = makePDB(
            ['Header'], atoms_list_4, ['Footer'], [], 'tests/temp_files/test.pdb', 'bfactor'
        )
        with open('tests/temp_files/test.pdb', 'r') as f:
            act_pdb_lines_4 = f.read().split('\n')
        exp_pdb_lines_4 = ['Header', 'Footer', '']
        self.assertEqual(act_pdb_lines_4, exp_pdb_lines_4)
        self.assertTrue(exit_4)

        atoms_list_5 = copy.deepcopy(atoms_list_1)
        atoms_list_5[0].resiType = 'GLUU'
        exit_5 = makePDB(
            ['Header'], atoms_list_5, ['Footer'], [], 'tests/temp_files/test.pdb', 'bfactor'
        )
        with open('tests/temp_files/test.pdb', 'r') as f:
            act_pdb_lines_5 = f.read().split('\n')
        exp_pdb_lines_5 = ['Header', 'Footer', '']
        self.assertEqual(act_pdb_lines_5, exp_pdb_lines_5)
        self.assertTrue(exit_5)

        atoms_list_6 = copy.deepcopy(atoms_list_1)
        atoms_list_6[0].chainID = 'AA'
        exit_6 = makePDB(
            ['Header'], atoms_list_6, ['Footer'], [], 'tests/temp_files/test.pdb', 'bfactor'
        )
        with open('tests/temp_files/test.pdb', 'r') as f:
            act_pdb_lines_6 = f.read().split('\n')
        exp_pdb_lines_6 = ['Header', 'Footer', '']
        self.assertEqual(act_pdb_lines_6, exp_pdb_lines_6)
        self.assertTrue(exit_6)

        atoms_list_7 = copy.deepcopy(atoms_list_1)
        atoms_list_7[0].insCode = 'AA'
        exit_7 = makePDB(
            ['Header'], atoms_list_7, ['Footer'], [], 'tests/temp_files/test.pdb', 'bfactor'
        )
        with open('tests/temp_files/test.pdb', 'r') as f:
            act_pdb_lines_7 = f.read().split('\n')
        exp_pdb_lines_7 = ['Header', 'Footer', '']
        self.assertEqual(act_pdb_lines_7, exp_pdb_lines_7)
        self.assertTrue(exit_7)

        atoms_list_8 = copy.deepcopy(atoms_list_1)
        atoms_list_8[0].xyzCoords[2][0] = 123456.7
        exit_8 = makePDB(
            ['Header'], atoms_list_8, ['Footer'], [], 'tests/temp_files/test.pdb', 'bfactor'
        )
        with open('tests/temp_files/test.pdb', 'r') as f:
            act_pdb_lines_8 = f.read().split('\n')
        exp_pdb_lines_8 = ['Header', 'Footer', '']
        self.assertEqual(act_pdb_lines_8, exp_pdb_lines_8)
        self.assertTrue(exit_8)

        atoms_list_9 = copy.deepcopy(atoms_list_1)
        atoms_list_9[0].occupancy = 12345
        exit_9 = makePDB(
            ['Header'], atoms_list_9, ['Footer'], [], 'tests/temp_files/test.pdb', 'bfactor'
        )
        with open('tests/temp_files/test.pdb', 'r') as f:
            act_pdb_lines_9 = f.read().split('\n')
        exp_pdb_lines_9 = ['Header', 'Footer', '']
        self.assertEqual(act_pdb_lines_9, exp_pdb_lines_9)
        self.assertTrue(exit_9)

        atoms_list_10 = copy.deepcopy(atoms_list_1)
        atoms_list_10[0].bFactor = 12345
        exit_10 = makePDB(
            ['Header'], atoms_list_10, ['Footer'], [], 'tests/temp_files/test.pdb', 'bfactor'
        )
        with open('tests/temp_files/test.pdb', 'r') as f:
            act_pdb_lines_10 = f.read().split('\n')
        exp_pdb_lines_10 = ['Header', 'Footer', '']
        self.assertEqual(act_pdb_lines_10, exp_pdb_lines_10)
        self.assertTrue(exit_10)

        atoms_list_11 = copy.deepcopy(atoms_list_1)
        atoms_list_11[0].element = 'FEE'
        exit_11 = makePDB(
            ['Header'], atoms_list_11, ['Footer'], [], 'tests/temp_files/test.pdb', 'bfactor'
        )
        with open('tests/temp_files/test.pdb', 'r') as f:
            act_pdb_lines_11 = f.read().split('\n')
        exp_pdb_lines_11 = ['Header', 'Footer', '']
        self.assertEqual(act_pdb_lines_11, exp_pdb_lines_11)
        self.assertTrue(exit_11)

        atoms_list_12 = copy.deepcopy(atoms_list_1)
        atoms_list_12[0].charge = '++2'
        exit_12 = makePDB(
            ['Header'], atoms_list_12, ['Footer'], [], 'tests/temp_files/test.pdb', 'bfactor'
        )
        with open('tests/temp_files/test.pdb', 'r') as f:
            act_pdb_lines_12 = f.read().split('\n')
        exp_pdb_lines_12 = ['Header', 'Footer', '']
        self.assertEqual(act_pdb_lines_12, exp_pdb_lines_12)
        self.assertTrue(exit_12)

        exit_13 = makePDB(
            ['Header'], atoms_list_1, ['Footer'], ['GLY', 'HOH', 'FAD'], 'tests/temp_files/test.pdb', 'bfactor'
        )
        with open('tests/temp_files/test.pdb', 'r') as f:
            act_pdb_lines_13 = f.read().split('\n')
        exp_pdb_lines_13 = [
            'Header',
            'ATOM      2  CA  GLY A   1      13.602  45.768  30.728  1.00 19.08           C  ',
            'TER                                                                             ',
            'HETATM  438  O   HOH B2001      13.128  46.298  34.342  1.00 33.10           O  ',
            'HETATM 3546  PA  FAD B 700      29.643  29.700  51.402  1.00 20.52           P  ',
            'TER                                                                             ',
            'Footer',
            ''
        ]
        self.assertEqual(act_pdb_lines_13, exp_pdb_lines_13)
        self.assertFalse(exit_13)

        exit_14 = makePDB(
            ['Header'], atoms_list_1, ['Footer'], [], 'tests/temp_files/test.pdb', 'bdamage'
        )
        with open('tests/temp_files/test.pdb', 'r') as f:
            act_pdb_lines_14 = f.read().split('\n')
        exp_pdb_lines_14 = [
            'Header',
            'ATOM      2  CA  GLY A   1      13.602  45.768  30.728  1.00  2.95           C  ',
            'HETATM  438  O   HOH B2001      13.128  46.298  34.342  1.00  3.50           O  ',
            'HETATM 3546  PA  FAD B 700      29.643  29.700  51.402  1.00  3.02           P  ',
            'Footer',
            ''
        ]
        self.assertEqual(act_pdb_lines_14, exp_pdb_lines_14)
        self.assertFalse(exit_14)

        shutil.rmtree('tests/temp_files/')
