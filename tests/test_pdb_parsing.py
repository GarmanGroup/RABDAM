
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


# python -m unittest tests/test_pdb_parsing.py

import unittest
from rabdam.Subroutines.PDBCUR import (
    get_rfree_from_pdb, get_res_temp_from_pdb, parse_atom_rec_from_pdb
)
from rabdam.Subroutines.makeDataFrame import makePDB, make_c_pdb
from tests.make_atom_rec import (
    gen_atom_objs_list, make_exp_pdb_rec
)


class TestClass(unittest.TestCase):

    def test_get_rfree_from_pdb(self):
        """
        Checks that Rfree value is parsed correctly from an input PDB file.
        Examples are taken from
        https://www.wwpdb.org/documentation/file-format-content/format33/remark3.html
        """

        # Check parses remarks from REFMAC correctly
        remark_rec_1 = [
            'REMARK   3  FIT TO DATA USED IN REFINEMENT.',
            'REMARK   3   CROSS-VALIDATION METHOD          : THROUGHOUT',
            'REMARK   3   FREE R VALUE TEST SET SELECTION  : RANDOM',
            'REMARK   3   R VALUE     (WORKING + TEST SET) : 0.228',
            'REMARK   3   R VALUE            (WORKING SET) : 0.225',
            'REMARK   3   FREE R VALUE                     : 0.283',
            'REMARK   3   FREE R VALUE TEST SET SIZE   (%) : 5.200',
            'REMARK   3   FREE R VALUE TEST SET COUNT      : 2256',
            'REMARK   3',
            'REMARK   3  FIT IN THE HIGHEST RESOLUTION BIN.',
            'REMARK   3   TOTAL NUMBER OF BINS USED           : 20',
            'REMARK   3   BIN RESOLUTION RANGE HIGH       (A) : 2.20',
            'REMARK   3   BIN RESOLUTION RANGE LOW        (A) : 2.26',
            'REMARK   3   REFLECTION IN BIN     (WORKING SET) : 2978',
            'REMARK   3   BIN COMPLETENESS (WORKING+TEST) (%) : 98.65',
            'REMARK   3   BIN R VALUE           (WORKING SET) : 0.2840',
            'REMARK   3   BIN FREE R VALUE SET COUNT          : 161',
            'REMARK   3   BIN FREE R VALUE                    : 0.3680'
        ]
        act_rfree_1 = get_rfree_from_pdb(remark_rec_1)
        exp_rfree_1 = 0.283
        self.assertEqual(round(act_rfree_1, 3), round(exp_rfree_1, 3))


        # Check parses remarks from SHELX correctly
        remark_rec_2 = [
            'REMARK   3  DATA USED IN REFINEMENT.',
            'REMARK   3   RESOLUTION RANGE HIGH (ANGSTROMS) : 1.15',
            'REMARK   3   RESOLUTION RANGE LOW  (ANGSTROMS) : 30.00',
            'REMARK   3   DATA CUTOFF            (SIGMA(F)) : 0.000',
            'REMARK   3   COMPLETENESS FOR RANGE        (%) : 99.8',
            'REMARK   3   CROSS-VALIDATION METHOD           : FREE R',
            'REMARK   3   FREE R VALUE TEST SET SELECTION   : RANDOM',
            'REMARK   3',
            'REMARK   3  FIT TO DATA USED IN REFINEMENT (NO CUTOFF).',
            'REMARK   3   R VALUE   (WORKING + TEST SET, NO CUTOFF) : 0.116',
            'REMARK   3   R VALUE          (WORKING SET, NO CUTOFF) : 0.116',
            'REMARK   3   FREE R VALUE                  (NO CUTOFF) : 0.145',
            'REMARK   3   FREE R VALUE TEST SET SIZE (%, NO CUTOFF) : 5.000',
            'REMARK   3   FREE R VALUE TEST SET COUNT   (NO CUTOFF) : 4279',
            'REMARK   3   TOTAL NUMBER OF REFLECTIONS   (NO CUTOFF) : 85756',
            'REMARK   3',
            'REMARK   3  FIT/AGREEMENT OF MODEL FOR DATA WITH F>4SIG(F).',
            'REMARK   3   R VALUE   (WORKING + TEST SET, F>4SIG(F)) : 0.010',
            'REMARK   3   R VALUE          (WORKING SET, F>4SIG(F)) : 0.010',
            'REMARK   3   FREE R VALUE                  (F>4SIG(F)) : 0.136',
            'REMARK   3   FREE R VALUE TEST SET SIZE (%, F>4SIG(F)) : 5.000',
            'REMARK   3   FREE R VALUE TEST SET COUNT   (F>4SIG(F)) : 3859',
            'REMARK   3   TOTAL NUMBER OF REFLECTIONS   (F>4SIG(F)) : 77074'
        ]
        act_rfree_2 = get_rfree_from_pdb(remark_rec_2)
        exp_rfree_2 = 0.145
        self.assertEqual(round(act_rfree_2, 3), round(exp_rfree_2, 3))

        # Check parses remarks from Phenix correctly
        remark_rec_3 = [
            'REMARK   3  DATA USED IN REFINEMENT.                                            ',
            'REMARK   3   RESOLUTION RANGE HIGH (ANGSTROMS) : 2.99                           ',
            'REMARK   3   RESOLUTION RANGE LOW  (ANGSTROMS) : 40.07                          ',
            'REMARK   3   MIN(FOBS/SIGMA_FOBS)              : 0.000                          ',
            'REMARK   3   COMPLETENESS FOR RANGE        (%) : 96.7                           ',
            'REMARK   3   NUMBER OF REFLECTIONS             : 242645                         ',
            'REMARK   3                                                                      ',
            'REMARK   3  FIT TO DATA USED IN REFINEMENT.                                     ',
            'REMARK   3   R VALUE     (WORKING + TEST SET) : 0.293                           ',
            'REMARK   3   R VALUE            (WORKING SET) : 0.291                           ',
            'REMARK   3   FREE R VALUE                     : 0.335                           ',
            'REMARK   3   FREE R VALUE TEST SET SIZE   (%) : 4.980                           ',
            'REMARK   3   FREE R VALUE TEST SET COUNT      : 12081                           '
        ]
        act_rfree_3 = get_rfree_from_pdb(remark_rec_3)
        exp_rfree_3 = 0.335
        self.assertEqual(round(act_rfree_3, 3), round(exp_rfree_3, 3))

        # Check returns None if Rfree isn't recorded
        remark_rec_4 = [
            'REMARK   3  DATA USED IN REFINEMENT.                                            ',
            'REMARK   3   RESOLUTION RANGE HIGH (ANGSTROMS) : 2.99                           ',
            'REMARK   3   RESOLUTION RANGE LOW  (ANGSTROMS) : 40.07                          ',
            'REMARK   3   MIN(FOBS/SIGMA_FOBS)              : 0.000                          ',
            'REMARK   3   COMPLETENESS FOR RANGE        (%) : 96.7                           ',
            'REMARK   3   NUMBER OF REFLECTIONS             : 242645                         ',
            'REMARK   3                                                                      ',
            'REMARK   3  FIT TO DATA USED IN REFINEMENT.                                     ',
            'REMARK   3   R VALUE     (WORKING + TEST SET) : 0.293                           ',
            'REMARK   3   R VALUE            (WORKING SET) : 0.291                           ',
            'REMARK   3   FREE R VALUE                     :                                 ',
            'REMARK   3   FREE R VALUE TEST SET SIZE   (%) : 4.980                           ',
            'REMARK   3   FREE R VALUE TEST SET COUNT      : 12081                           '
        ]
        act_rfree_4 = get_rfree_from_pdb(remark_rec_4)
        self.assertIsNone(act_rfree_4)

        # Check returns None if Rfree remark isn't included
        remark_rec_5 = [
            'REMARK   3  DATA USED IN REFINEMENT.                                            ',
            'REMARK   3   RESOLUTION RANGE HIGH (ANGSTROMS) : 2.99                           ',
            'REMARK   3   RESOLUTION RANGE LOW  (ANGSTROMS) : 40.07                          ',
            'REMARK   3   MIN(FOBS/SIGMA_FOBS)              : 0.000                          ',
            'REMARK   3   COMPLETENESS FOR RANGE        (%) : 96.7                           ',
            'REMARK   3   NUMBER OF REFLECTIONS             : 242645                         ',
            'REMARK   3                                                                      ',
            'REMARK   3  FIT TO DATA USED IN REFINEMENT.                                     ',
            'REMARK   3   R VALUE     (WORKING + TEST SET) : 0.293                           ',
            'REMARK   3   R VALUE            (WORKING SET) : 0.291                           ',
            'REMARK   3   FREE R VALUE TEST SET SIZE   (%) : 4.980                           ',
            'REMARK   3   FREE R VALUE TEST SET COUNT      : 12081                           '
        ]
        act_rfree_5 = get_rfree_from_pdb(remark_rec_5)
        self.assertIsNone(act_rfree_5)

        # Check returns None if Rfree remark is set to a non-numeric value
        remark_rec_6 = [
            'REMARK   3  DATA USED IN REFINEMENT.                                            ',
            'REMARK   3   RESOLUTION RANGE HIGH (ANGSTROMS) : 2.99                           ',
            'REMARK   3   RESOLUTION RANGE LOW  (ANGSTROMS) : 40.07                          ',
            'REMARK   3   MIN(FOBS/SIGMA_FOBS)              : 0.000                          ',
            'REMARK   3   COMPLETENESS FOR RANGE        (%) : 96.7                           ',
            'REMARK   3   NUMBER OF REFLECTIONS             : 242645                         ',
            'REMARK   3                                                                      ',
            'REMARK   3  FIT TO DATA USED IN REFINEMENT.                                     ',
            'REMARK   3   R VALUE     (WORKING + TEST SET) : 0.293                           ',
            'REMARK   3   R VALUE            (WORKING SET) : 0.291                           ',
            'REMARK   3   FREE R VALUE                     : N/A                             ',
            'REMARK   3   FREE R VALUE TEST SET SIZE   (%) : 4.980                           ',
            'REMARK   3   FREE R VALUE TEST SET COUNT      : 12081                           '
        ]
        act_rfree_6 = get_rfree_from_pdb(remark_rec_6)
        self.assertIsNone(act_rfree_6)

    def test_get_res_temp_from_pdb(self):
        """
        Checks that resolution and temperature values are parsed correctly from
        an input PDB file.
        REMARK 3 records are taken from
        Examples are taken from
        https://www.wwpdb.org/documentation/file-format-content/format33/remark3.html.
        REMARK 200 records are taken from
        https://www.wwpdb.org/documentation/file-format-content/format33/remarks1.html.
        """

        # Check resolution is parsed correctly
        remark_rec_1 = [
            'REMARK   3  DATA USED IN REFINEMENT.',
            'REMARK   3   RESOLUTION RANGE HIGH (ANGSTROMS) : 1.15',
            'REMARK   3   RESOLUTION RANGE LOW  (ANGSTROMS) : 30.00',
            'REMARK   3   DATA CUTOFF            (SIGMA(F)) : 0.000',
            'REMARK   3   COMPLETENESS FOR RANGE        (%) : 99.8',
            'REMARK   3   CROSS-VALIDATION METHOD           : FREE R',
            'REMARK   3   FREE R VALUE TEST SET SELECTION   : RANDOM',
        ]
        act_resolution_1, act_temperature_1 = get_res_temp_from_pdb(remark_rec_1)
        exp_resolution_1 = 1.15
        self.assertEqual(round(act_resolution_1, 2), round(exp_resolution_1, 2))
        self.assertIsNone(act_temperature_1)

        # Check temperature is parsed correctly
        remark_rec_2 = [
            'REMARK 200                                                                      ',
            'REMARK 200 EXPERIMENTAL DETAILS                                                 ',
            'REMARK 200  EXPERIMENT TYPE                : X-RAY DIFFRACTION                  ',
            'REMARK 200  DATE OF  DATA COLLECTION       : 17-MAR-02;  17-MAR-02              ',
            'REMARK 200  TEMPERATURE           (KELVIN) : 100                                ',
            'REMARK 200  PH                             : 8.00                               ',
            'REMARK 200  NUMBER  OF CRYSTALS USED       : 2                                  ',
        ]
        act_resolution_2, act_temperature_2 = get_res_temp_from_pdb(remark_rec_2)
        exp_temperature_2 = [100.0]
        self.assertEqual(act_temperature_2, exp_temperature_2)
        self.assertIsNone(act_resolution_2)

        # Check temperature list is parsed correctly
        remark_rec_3 = [
            'REMARK 200                                                                      ',
            'REMARK 200 EXPERIMENTAL DETAILS                                                 ',
            'REMARK 200  EXPERIMENT TYPE                : X-RAY DIFFRACTION                  ',
            'REMARK 200  DATE OF  DATA COLLECTION       : 17-MAR-02;  17-MAR-02              ',
            'REMARK 200  TEMPERATURE           (KELVIN) : 100; 100                           ',
            'REMARK 200  PH                             : 8.00                               ',
            'REMARK 200  NUMBER  OF CRYSTALS USED       : 2                                  ',
        ]
        act_resolution_3, act_temperature_3 = get_res_temp_from_pdb(remark_rec_3)
        exp_temperature_3 = [100.0, 100.0]
        self.assertEqual(act_temperature_3, exp_temperature_3)
        self.assertIsNone(act_resolution_3)

        # Check non-numeric resolution is parsed correctly
        remark_rec_4 = [
            'REMARK   3  DATA USED IN REFINEMENT.',
            'REMARK   3   RESOLUTION RANGE HIGH (ANGSTROMS) : N/A',
            'REMARK   3   RESOLUTION RANGE LOW  (ANGSTROMS) : 30.00',
            'REMARK   3   DATA CUTOFF            (SIGMA(F)) : 0.000',
            'REMARK   3   COMPLETENESS FOR RANGE        (%) : 99.8',
            'REMARK   3   CROSS-VALIDATION METHOD           : FREE R',
            'REMARK   3   FREE R VALUE TEST SET SELECTION   : RANDOM',
        ]
        act_resolution_4, act_temperature_4 = get_res_temp_from_pdb(remark_rec_4)
        self.assertIsNone(act_temperature_4)
        self.assertIsNone(act_resolution_4)

        # Check non-numeric temperature is parsed correctly
        remark_rec_5 = [
            'REMARK 200                                                                      ',
            'REMARK 200 EXPERIMENTAL DETAILS                                                 ',
            'REMARK 200  EXPERIMENT TYPE                : X-RAY DIFFRACTION                  ',
            'REMARK 200  DATE OF  DATA COLLECTION       : 17-MAR-02;  17-MAR-02              ',
            'REMARK 200  TEMPERATURE           (KELVIN) : 100; N/A                           ',
            'REMARK 200  PH                             : 8.00                               ',
            'REMARK 200  NUMBER  OF CRYSTALS USED       : 2                                  ',
        ]
        act_resolution_5, act_temperature_5 = get_res_temp_from_pdb(remark_rec_5)
        self.assertIsNone(act_temperature_5)
        self.assertIsNone(act_resolution_5)

        # Check unspecified resolution is parsed correctly
        remark_rec_6 = [
            'REMARK   3  DATA USED IN REFINEMENT.',
            'REMARK   3   RESOLUTION RANGE HIGH (ANGSTROMS) : ',
            'REMARK   3   RESOLUTION RANGE LOW  (ANGSTROMS) : 30.00',
            'REMARK   3   DATA CUTOFF            (SIGMA(F)) : 0.000',
            'REMARK   3   COMPLETENESS FOR RANGE        (%) : 99.8',
            'REMARK   3   CROSS-VALIDATION METHOD           : FREE R',
            'REMARK   3   FREE R VALUE TEST SET SELECTION   : RANDOM',
        ]
        act_resolution_6, act_temperature_6 = get_res_temp_from_pdb(remark_rec_6)
        self.assertIsNone(act_temperature_6)
        self.assertIsNone(act_resolution_6)

        # Check unspecified temperature is parsed correctly
        remark_rec_7 = [
            'REMARK 200                                                                      ',
            'REMARK 200 EXPERIMENTAL DETAILS                                                 ',
            'REMARK 200  EXPERIMENT TYPE                : X-RAY DIFFRACTION                  ',
            'REMARK 200  DATE OF  DATA COLLECTION       : 17-MAR-02;  17-MAR-02              ',
            'REMARK 200  TEMPERATURE           (KELVIN) :                                    ',
            'REMARK 200  PH                             : 8.00                               ',
            'REMARK 200  NUMBER  OF CRYSTALS USED       : 2                                  ',
        ]
        act_resolution_7, act_temperature_7 = get_res_temp_from_pdb(remark_rec_7)
        self.assertIsNone(act_temperature_7)
        self.assertIsNone(act_resolution_7)

    def test_parse_atomrec_from_pdb(self):
        """
        Checks that PDB format files are parsed as expected
        """

        import copy
        import os
        import shutil

        # Load example PDB records
        exp_atom_rec_obj_1, _, exp_atom_rec_str_1 = make_exp_pdb_rec()
        if not os.path.isdir('tests/temp_files/'):
            os.mkdir('tests/temp_files/')
        with open('tests/temp_files/test.pdb', 'w') as f:
            f.write('\n'.join(exp_atom_rec_str_1))

        # Check that atom records are parsed correctly
        act_atom_rec_1, exit_1 = parse_atom_rec_from_pdb(
            exp_atom_rec_str_1, 'tests/temp_files/test.pdb', False, False
        )
        for index in range(len(act_atom_rec_1)):
            self.assertEqual(exp_atom_rec_obj_1[index], act_atom_rec_1[index])
        self.assertFalse(exit_1)

        # Check raises error if resiNum is not a numeric value
        exp_atom_rec_str_2 = copy.deepcopy(exp_atom_rec_str_1)
        exp_atom_rec_str_2[0] = (
            exp_atom_rec_str_2[0][:22] + '   ?' + exp_atom_rec_str_2[0][26:]
        )
        act_atom_rec_2, exit_2 = parse_atom_rec_from_pdb(
            exp_atom_rec_str_2, 'tests/temp_files/test.pdb', False, False
        )
        self.assertTrue(exit_2)

        # Check raises error if atomNum is not a numeric value
        exp_atom_rec_str_3 = copy.deepcopy(exp_atom_rec_str_1)
        exp_atom_rec_str_3[0] = (
            exp_atom_rec_str_3[0][:6] + '    ?' + exp_atom_rec_str_3[0][11:]
        )
        act_atom_rec_3, exit_3 = parse_atom_rec_from_pdb(
            exp_atom_rec_str_3, 'tests/temp_files/test.pdb', False, False
        )
        self.assertTrue(exit_3)

        # Check raises error if xyzCoords are not numeric values
        exp_atom_rec_str_4 = copy.deepcopy(exp_atom_rec_str_1)
        exp_atom_rec_str_4[0] = (
            exp_atom_rec_str_4[0][:30] + '       ?' + exp_atom_rec_str_4[0][38:]
        )
        act_atom_rec_4, exit_4 = parse_atom_rec_from_pdb(
            exp_atom_rec_str_4, 'tests/temp_files/test.pdb', False, False
        )
        self.assertTrue(exit_4)

        # Check raises error if occupancy is not a numeric value
        exp_atom_rec_str_5 = copy.deepcopy(exp_atom_rec_str_1)
        exp_atom_rec_str_5[0] = (
            exp_atom_rec_str_5[0][:54] + '     ?' + exp_atom_rec_str_5[0][60:]
        )
        act_atom_rec_5, exit_5 = parse_atom_rec_from_pdb(
            exp_atom_rec_str_5, 'tests/temp_files/test.pdb', False, False
        )
        self.assertTrue(exit_5)

        # Check raises error if bFactor is not a numeric value
        exp_atom_rec_str_6 = copy.deepcopy(exp_atom_rec_str_1)
        exp_atom_rec_str_6[0] = (
            exp_atom_rec_str_6[0][:60] + '     ?' + exp_atom_rec_str_6[0][66:]
        )
        act_atom_rec_6, exit_6 = parse_atom_rec_from_pdb(
            exp_atom_rec_str_6, 'tests/temp_files/test.pdb', False, False
        )
        self.assertTrue(exit_6)

        # Check raises error if record columns are missing
        exp_atom_rec_str_7 = [
            line[:76] for line in copy.deepcopy(exp_atom_rec_str_1)
        ]
        act_atom_rec_7, exit_7 = parse_atom_rec_from_pdb(
            exp_atom_rec_str_7, 'tests/temp_files/test.pdb', False, False
        )
        self.assertTrue(exit_7)

    def test_write_pdb_file(self):
        """
        Tests that atom objects are correctly converted into PDB format, and TER
        cards are inserted at appropriate locations
        """

        import numpy as np
        import os
        import copy
        import shutil

        if not os.path.isdir('tests/temp_files/'):
            os.mkdir('tests/temp_files/')

        # Generate list of atom objects
        atoms_list_1 = gen_atom_objs_list()

        # Check PDB file is created as expected with B-factors
        exit_1 = makePDB(
            ['Header'], atoms_list_1, ['Footer'], 'tests/temp_files/test.pdb',
            'bfactor'
        )
        with open('tests/temp_files/test.pdb', 'r') as f:
            act_pdb_lines_1 = f.read().split('\n')
        exp_pdb_lines_1 = [
            'Header',
            'ATOM      2  CA  GLY A   1      13.602  45.768  30.728  1.00 19.08           C  ',
            'HETATM  438  O   HOH B2001      13.128  46.298  34.342  1.00 33.10           O  ',
            'HETATM 3546  PA XFAD B 700      29.643  29.700  51.402  1.00 20.52           P  ',
            'Footer',
            'END                                                                             ',
            ''
        ]
        self.assertEqual(act_pdb_lines_1, exp_pdb_lines_1)
        self.assertFalse(exit_1)

        # Check PDB file is created as expected with BDamage values in the
        # B-factor column
        exit_2 = makePDB(
            ['Header'], atoms_list_1, ['Footer'], 'tests/temp_files/test.pdb',
            'bdamage'
        )
        with open('tests/temp_files/test.pdb', 'r') as f:
            act_pdb_lines_2 = f.read().split('\n')
        exp_pdb_lines_2 = [
            'Header',
            'ATOM      2  CA  GLY A   1      13.602  45.768  30.728  1.00  0.26           C  ',
            'HETATM  438  O   HOH B2001      13.128  46.298  34.342  1.00 -0.03           O  ',
            'HETATM 3546  PA XFAD B 700      29.643  29.700  51.402  1.00  0.69           P  ',
            'Footer',
            'END                                                                             ',
            ''
        ]
        self.assertEqual(act_pdb_lines_2, exp_pdb_lines_2)
        self.assertFalse(exit_2)

        # Check raises error if lineID is longer than 6 characters
        atoms_list_3 = copy.deepcopy(atoms_list_1)
        atoms_list_3[0].lineID = 'HETATMM'
        exit_3 = makePDB(
            ['Header'], atoms_list_3, ['Footer'], 'tests/temp_files/test.pdb',
            'bfactor'
        )
        with open('tests/temp_files/test.pdb', 'r') as f:
            act_pdb_lines_3 = f.read().split('\n')
        exp_pdb_lines_3 = [
            'Header',
            'Footer',
            'END                                                                             ',
            ''
        ]
        self.assertEqual(act_pdb_lines_3, exp_pdb_lines_3)
        self.assertTrue(exit_3)

        # Check raises error if atomType is longer than 3 characters
        atoms_list_4 = copy.deepcopy(atoms_list_1)
        atoms_list_4[0].atomType = 'HA12'
        exit_4 = makePDB(
            ['Header'], atoms_list_4, ['Footer'], 'tests/temp_files/test.pdb',
            'bfactor'
        )
        with open('tests/temp_files/test.pdb', 'r') as f:
            act_pdb_lines_4 = f.read().split('\n')
        exp_pdb_lines_4 = [
            'Header',
            'Footer',
            'END                                                                             ',
            ''
        ]
        self.assertEqual(act_pdb_lines_4, exp_pdb_lines_4)
        self.assertTrue(exit_4)

        # Check raises error if conformer is longer than 1 character
        atoms_list_5 = copy.deepcopy(atoms_list_1)
        atoms_list_5[0].conformer = 'AA'
        exit_5 = makePDB(
            ['Header'], atoms_list_5, ['Footer'], 'tests/temp_files/test.pdb',
            'bfactor'
        )
        with open('tests/temp_files/test.pdb', 'r') as f:
            act_pdb_lines_5 = f.read().split('\n')
        exp_pdb_lines_5 = [
            'Header',
            'Footer',
            'END                                                                             ',
            ''
        ]
        self.assertEqual(act_pdb_lines_5, exp_pdb_lines_5)
        self.assertTrue(exit_5)

        # Check raises error if resiType is longer than 3 characters
        atoms_list_6 = copy.deepcopy(atoms_list_1)
        atoms_list_6[0].resiType = 'GLYY'
        exit_6 = makePDB(
            ['Header'], atoms_list_6, ['Footer'], 'tests/temp_files/test.pdb',
            'bfactor'
        )
        with open('tests/temp_files/test.pdb', 'r') as f:
            act_pdb_lines_6 = f.read().split('\n')
        exp_pdb_lines_6 = [
            'Header',
            'Footer',
            'END                                                                             ',
            ''
        ]
        self.assertEqual(act_pdb_lines_6, exp_pdb_lines_6)
        self.assertTrue(exit_6)

        # Check raises error if chainID is longer than 1 character
        atoms_list_7 = copy.deepcopy(atoms_list_1)
        atoms_list_7[0].chainID = 'AA'
        exit_7 = makePDB(
            ['Header'], atoms_list_7, ['Footer'], 'tests/temp_files/test.pdb',
            'bfactor'
        )
        with open('tests/temp_files/test.pdb', 'r') as f:
            act_pdb_lines_7 = f.read().split('\n')
        exp_pdb_lines_7 = [
            'Header',
            'Footer',
            'END                                                                             ',
            ''
        ]
        self.assertEqual(act_pdb_lines_7, exp_pdb_lines_7)
        self.assertTrue(exit_7)

        # Check raises error if insCode is longer than 1 character
        atoms_list_8 = copy.deepcopy(atoms_list_1)
        atoms_list_8[0].insCode = 'AA'
        exit_8 = makePDB(
            ['Header'], atoms_list_8, ['Footer'], 'tests/temp_files/test.pdb',
            'bfactor'
        )
        with open('tests/temp_files/test.pdb', 'r') as f:
            act_pdb_lines_8 = f.read().split('\n')
        exp_pdb_lines_8 = [
            'Header',
            'Footer',
            'END                                                                             ',
            ''
        ]
        self.assertEqual(act_pdb_lines_8, exp_pdb_lines_8)
        self.assertTrue(exit_8)

        # Check rounds xyzCoords correctly
        atoms_list_9 = copy.deepcopy(atoms_list_1)
        atoms_list_9[0].xyzCoords = np.array([[1234.5678], [1234.5678], [1234.5678]])
        exit_9 = makePDB(
            ['Header'], atoms_list_9, ['Footer'], 'tests/temp_files/test.pdb',
            'bdamage'
        )
        with open('tests/temp_files/test.pdb', 'r') as f:
            act_pdb_lines_9 = f.read().split('\n')
        exp_pdb_lines_9 = [
            'Header',
            'ATOM      2  CA  GLY A   1    1234.5681234.5681234.568  1.00  0.26           C  ',
            'HETATM  438  O   HOH B2001      13.128  46.298  34.342  1.00 -0.03           O  ',
            'HETATM 3546  PA XFAD B 700      29.643  29.700  51.402  1.00  0.69           P  ',
            'Footer',
            'END                                                                             ',
            ''
        ]
        self.assertEqual(act_pdb_lines_9, exp_pdb_lines_9)
        self.assertFalse(exit_9)

        # Check raises error if xyzCoords are longer than 8 characters
        atoms_list_10 = copy.deepcopy(atoms_list_1)
        atoms_list_10[0].xyzCoords = np.array([[12345.5678], [1234.5678], [1234.5678]])
        exit_10 = makePDB(
            ['Header'], atoms_list_10, ['Footer'], 'tests/temp_files/test.pdb',
            'bfactor'
        )
        with open('tests/temp_files/test.pdb', 'r') as f:
            act_pdb_lines_10 = f.read().split('\n')
        exp_pdb_lines_10 = [
            'Header',
            'Footer',
            'END                                                                             ',
            ''
        ]
        self.assertEqual(act_pdb_lines_10, exp_pdb_lines_10)
        self.assertTrue(exit_10)

        # Check rounds occupancy correctly
        atoms_list_11 = copy.deepcopy(atoms_list_1)
        atoms_list_11[0].occupancy = 1.2345
        exit_11 = makePDB(
            ['Header'], atoms_list_11, ['Footer'], 'tests/temp_files/test.pdb',
            'bdamage'
        )
        with open('tests/temp_files/test.pdb', 'r') as f:
            act_pdb_lines_11 = f.read().split('\n')
        exp_pdb_lines_11 = [
            'Header',
            'ATOM      2  CA  GLY A   1      13.602  45.768  30.728  1.23  0.26           C  ',
            'HETATM  438  O   HOH B2001      13.128  46.298  34.342  1.00 -0.03           O  ',
            'HETATM 3546  PA XFAD B 700      29.643  29.700  51.402  1.00  0.69           P  ',
            'Footer',
            'END                                                                             ',
            ''
        ]
        self.assertEqual(act_pdb_lines_11, exp_pdb_lines_11)
        self.assertFalse(exit_11)

        # Check raises error if occupancy is longer than 6 characters
        atoms_list_12 = copy.deepcopy(atoms_list_1)
        atoms_list_12[0].occupancy = 1234.56
        exit_12 = makePDB(
            ['Header'], atoms_list_12, ['Footer'], 'tests/temp_files/test.pdb',
            'bfactor'
        )
        with open('tests/temp_files/test.pdb', 'r') as f:
            act_pdb_lines_12 = f.read().split('\n')
        exp_pdb_lines_12 = [
            'Header',
            'Footer',
            'END                                                                             ',
            ''
        ]
        self.assertEqual(act_pdb_lines_12, exp_pdb_lines_12)
        self.assertTrue(exit_12)

        # Check rounds bFactor correctly
        atoms_list_13 = copy.deepcopy(atoms_list_1)
        atoms_list_13[0].bFactor = 123.456
        exit_13 = makePDB(
            ['Header'], atoms_list_13, ['Footer'], 'tests/temp_files/test.pdb',
            'bfactor'
        )
        with open('tests/temp_files/test.pdb', 'r') as f:
            act_pdb_lines_13 = f.read().split('\n')
        exp_pdb_lines_13 = [
            'Header',
            'ATOM      2  CA  GLY A   1      13.602  45.768  30.728  1.00123.46           C  ',
            'HETATM  438  O   HOH B2001      13.128  46.298  34.342  1.00 33.10           O  ',
            'HETATM 3546  PA XFAD B 700      29.643  29.700  51.402  1.00 20.52           P  ',
            'Footer',
            'END                                                                             ',
            ''
        ]
        self.assertEqual(act_pdb_lines_13, exp_pdb_lines_13)
        self.assertFalse(exit_13)

        # Check raises error if occupancy is longer than 6 characters
        atoms_list_14 = copy.deepcopy(atoms_list_1)
        atoms_list_14[0].bFactor = 1234.56
        exit_14 = makePDB(
            ['Header'], atoms_list_14, ['Footer'], 'tests/temp_files/test.pdb',
            'bfactor'
        )
        with open('tests/temp_files/test.pdb', 'r') as f:
            act_pdb_lines_14 = f.read().split('\n')
        exp_pdb_lines_14 = [
            'Header',
            'Footer',
            'END                                                                             ',
            ''
        ]
        self.assertEqual(act_pdb_lines_14, exp_pdb_lines_14)
        self.assertTrue(exit_14)

        # Check raises error if Bfac variable is set to something other than
        # 'bfactor' or 'bdamage'
        atoms_list_15 = copy.deepcopy(atoms_list_1)
        self.assertRaises(
            ValueError, makePDB, ['Header'], atoms_list_15, ['Footer'],
            'tests/temp_files/test.pdb', 'bfac'
        )

        # Check raises error if element is longer than 2 characters
        atoms_list_16 = copy.deepcopy(atoms_list_1)
        atoms_list_16[0].element = 'AAA'
        exit_16 = makePDB(
            ['Header'], atoms_list_16, ['Footer'], 'tests/temp_files/test.pdb',
            'bfactor'
        )
        with open('tests/temp_files/test.pdb', 'r') as f:
            act_pdb_lines_16 = f.read().split('\n')
        exp_pdb_lines_16 = [
            'Header',
            'Footer',
            'END                                                                             ',
            ''
        ]
        self.assertEqual(act_pdb_lines_16, exp_pdb_lines_16)
        self.assertTrue(exit_16)

        # Check raises error if charge is longer than 2 characters
        atoms_list_17 = copy.deepcopy(atoms_list_1)
        atoms_list_17[0].charge = '+27'
        exit_17 = makePDB(
            ['Header'], atoms_list_17, ['Footer'], 'tests/temp_files/test.pdb',
            'bfactor'
        )
        with open('tests/temp_files/test.pdb', 'r') as f:
            act_pdb_lines_17 = f.read().split('\n')
        exp_pdb_lines_17 = [
            'Header',
            'Footer',
            'END                                                                             ',
            ''
        ]
        self.assertEqual(act_pdb_lines_17, exp_pdb_lines_17)
        self.assertTrue(exit_17)

        # Check inserts TER cards correctly
        exp_atom_rec_obj_18, exp_atom_rec_str_18, _ = make_exp_pdb_rec()
        exit_18 = makePDB(
            [], exp_atom_rec_obj_18, [], 'tests/temp_files/test.pdb',
            'bfactor'
        )
        with open('tests/temp_files/test.pdb', 'r') as f:
            act_pdb_lines_18 = f.read().split('\n')
        # Filter to remove MG records - RABDAM makes mistakes with some metal
        # ions, doesn't affect the running of the PDB file in e.g. PyMol
        act_pdb_lines_18 = [
            line for line in act_pdb_lines_18 if not line[76:78] == 'MG'
        ]
        exp_atom_rec_str_18 = [
            line for line in exp_atom_rec_str_18 if not line[76:78] == 'MG'
        ] + ['']
        self.assertEqual(act_pdb_lines_18, exp_atom_rec_str_18)
        self.assertFalse(exit_18)

        shutil.rmtree('tests/temp_files/')

    def test_write_c_pdb_file(self):
        """
        Tests that atom objects are correctly converted into PDB format, and TER
        cards are inserted at appropriate locations
        """

        import numpy as np
        import os
        import copy
        import shutil

        if not os.path.isdir('tests/temp_files/'):
            os.mkdir('tests/temp_files/')

        # Generate list of atom objects
        atoms_list_1 = np.array([[1.0, 2.4, 2.5],
                                 [6.334567, 1234.567, 1.7],
                                 [1.48, 26.2, 9]])
        atom_id_list = [1, 2, 51]

        # Check PDB file is created as expected
        exit_1 = make_c_pdb(
            ['Header'], atoms_list_1, atom_id_list, ['Footer'],
            'tests/temp_files/test.pdb',
        )
        with open('tests/temp_files/test.pdb', 'r') as f:
            act_pdb_lines_1 = f.read().split('\n')
        exp_pdb_lines_1 = [
            'Header',
            'ATOM      1  C   GLY A   1       1.000   2.400   2.500  1.00 50.00           C  ',
            'ATOM      2  C   GLY A   1       6.3351234.567   1.700  1.00 50.00           C  ',
            'ATOM     51  C   GLY A   1       1.480  26.200   9.000  1.00 50.00           C  ',
            'Footer',
            'END                                                                             ',
            ''
        ]
        self.assertEqual(act_pdb_lines_1, exp_pdb_lines_1)
        self.assertFalse(exit_1)

        # Check raises error if xyzCoords are longer than 8 characters
        atoms_list_2 = copy.deepcopy(atoms_list_1)
        atoms_list_2[0][0] = 12345.5678
        exit_2 = make_c_pdb(
            ['Header'], atoms_list_2, atom_id_list, ['Footer'],
            'tests/temp_files/test.pdb',
        )
        with open('tests/temp_files/test.pdb', 'r') as f:
            act_pdb_lines_2 = f.read().split('\n')
        exp_pdb_lines_2 = [
            'Header',
            'Footer',
            'END                                                                             ',
            ''
        ]
        self.assertEqual(act_pdb_lines_2, exp_pdb_lines_2)
        self.assertTrue(exit_2)

        shutil.rmtree('tests/temp_files/')
