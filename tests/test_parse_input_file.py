
# RABDAM
# Copyright (C) 2025 Garman Group, University of Oxford

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


# python -m unittest tests/test_parse_input_file.py

import unittest
from rabdam.rabdam import (
    ArgumentError, parse_command_line_arguments, parse_input_file_arguments
)

class TestClass(unittest.TestCase):

    def test_parse_command_line_arguments(self):
        """
        Checks that command line arguments are parsed as expected
        """

        import argparse

        # Check raises SystemExist if no arguments provided
        command_line_1 = []
        self.assertRaises(
            SystemExit, parse_command_line_arguments, command_line_1, True
        )

        # Check parses -i flag correctly
        command_line_2 = ['-i', 'tests/test_files/example input.txt']
        act_parsed_2 = parse_command_line_arguments(command_line_2, True)
        exp_parsed_2 = argparse.Namespace(
            dependencies=False, version=False,
            input='tests/test_files/example input.txt', pdb_or_mmcif_file=None,
            run='full'
        )
        self.assertEqual(act_parsed_2, exp_parsed_2)

        # Check raises FileNotFoundError if input file specified by -i flag
        # doesn't exist
        command_line_3 = ['-i', 'does_not_exist.txt']
        self.assertRaises(
            FileNotFoundError, parse_command_line_arguments, command_line_3,
            True
        )

        # Check parses -f flag correctly with one PDB accession code
        command_line_4 = ['-f', '2bn1']
        act_parsed_4 = parse_command_line_arguments(command_line_4, True)
        exp_parsed_4 = argparse.Namespace(
            dependencies=False, version=False, input=None,
            pdb_or_mmcif_file=['2bn1'], run='full'
        )
        self.assertEqual(act_parsed_4, exp_parsed_4)

        # Check parses -f flag correctly with multiple PDB accession codes
        command_line_5 = ['-f', '5eeu', '5eev', '5eew']
        act_parsed_5 = parse_command_line_arguments(command_line_5, True)
        exp_parsed_5 = argparse.Namespace(
            dependencies=False, version=False, input=None,
            pdb_or_mmcif_file=['5eeu', '5eev', '5eew'], run='full'
        )
        self.assertEqual(act_parsed_5, exp_parsed_5)

        # Check parses -r flag correctly
        command_line_6 = ['-r', 'df', '-f', '2bn1']
        act_parsed_6 = parse_command_line_arguments(command_line_6, True)
        exp_parsed_6 = argparse.Namespace(
            dependencies=False, version=False, input=None,
            pdb_or_mmcif_file=['2bn1'], run='df'
        )
        self.assertEqual(act_parsed_6, exp_parsed_6)

        # Check parses --dependencies flag correctly
        command_line_7 = ['--dependencies']
        act_parsed_7 = parse_command_line_arguments(command_line_7, True)
        exp_parsed_7 = argparse.Namespace(
            dependencies=True, version=False, input=None,
            pdb_or_mmcif_file=None, run=None
        )
        self.assertEqual(act_parsed_7, exp_parsed_7)

        # Check parses --version flag correctly
        command_line_8 = ['--version']
        exp_version_8 = '115.0.6'
        act_parsed_8, act_version_8 = parse_command_line_arguments(
            command_line_8, True, exp_version_8
        )
        exp_parsed_8 = argparse.Namespace(
            dependencies=False, version=True, input=None,
            pdb_or_mmcif_file=None, run=None
        )
        self.assertEqual(act_parsed_8, exp_parsed_8)
        self.assertEqual(act_version_8, exp_version_8)

    def test_parse_input_file(self):
        """
        Checks that program arguments listed in input file are parsed as expected
        """

        import os

        cwd = os.getcwd()

        # Check the default input variables returned are correct
        input_1 = []
        act_output_1 = parse_input_file_arguments(input_1)
        exp_output_1 = {'outputDir': cwd,
                        'batchRun': False,
                        'overwrite': False,
                        'outFiles': 'all',
                        'filter': False,
                        'temperature': None,
                        'resolution': None,
                        'PDT': 7.0,
                        'windowSize': 0.02,
                        'HETATM': False,
                        'removeAtoms': [],
                        'addAtoms': [],
                        'highlightAtoms': [],
                        'createOrigpdb': False,
                        'createAUpdb': False,
                        'createUCpdb': False,
                        'createAUCpdb': False,
                        'createTApdb': False}
        self.assertDictEqual(act_output_1, exp_output_1)

        # Check that non-default, valid input variables are parsed correctly
        input_2 = [
            'outputdir=tests/', 'batchcontinue=y', 'overwrite=n',
            'outfiles=bnet', 'filter=n', 'temperature=120K', 'resolution=1.73',
            'pdt=4', 'windowsize=0.5', 'hetatm=keep', 'removeatoms=2;5;7-9;HOH',
            'addatoms=NA;FE;45-49;3;FAD', 'highlightatoms=1-6',
            'createorigpdb=YES', 'createaupdb=No', 'createucpdb=true',
            'createaucpdb=False', 'createtapdb=True'
        ]
        act_output_2 = parse_input_file_arguments(input_2)
        exp_output_2 = {'outputDir': 'tests/',
                        'batchRun': True,
                        'overwrite': False,
                        'outFiles': 'bnet',
                        'filter': False,
                        'temperature': 120,
                        'resolution': 1.73,
                        'PDT': 4.0,
                        'windowSize': 0.5,
                        'HETATM': True,
                        'removeAtoms': ['2', '5', '7', '8', '9', 'HOH'],
                        'addAtoms': ['NA', 'FE', '45', '46', '47', '48', '49', '3', 'FAD'],
                        'highlightAtoms': ['1', '2', '3', '4', '5', '6'],
                        'createOrigpdb': True,
                        'createAUpdb': False,
                        'createUCpdb': True,
                        'createAUCpdb': False,
                        'createTApdb': True}
        self.assertDictEqual(act_output_2, exp_output_2)

        # Check raises FileNotFoundError if output directory doesn't exist
        input_3 = ['outputdir=does_not_exist/']
        self.assertRaises(FileNotFoundError, parse_input_file_arguments, input_3)

        # Check raises ArgumentError if value specified for batchcontinue isn't
        # recognised
        input_4 = ['batchcontinue=xyz']
        self.assertRaises(ArgumentError, parse_input_file_arguments, input_4)

        # Check raises ArgumentError if value specified for overwrite isn't
        # recognised
        input_5 = ['overwrite=xyz']
        self.assertRaises(ArgumentError, parse_input_file_arguments, input_5)

        # Check raises ArgumentError if value specified for outfiles isn't
        # recognised
        input_6 = ['outfiles=xyz']
        self.assertRaises(ArgumentError, parse_input_file_arguments, input_6)

        # Check raises ArgumentError if value specified for filter isn't
        # recognised
        input_7 = ['filter=xyz']
        self.assertRaises(ArgumentError, parse_input_file_arguments, input_7)

        # Check raises ArgumentError if value specified for temperature isn't
        # recognised
        input_8 = ['temperature=xyz']
        self.assertRaises(ArgumentError, parse_input_file_arguments, input_8)

        # Check temperature='cryo' is converted to temperature=100
        input_9 = ['temperature=cryo']
        act_output_9 = parse_input_file_arguments(input_9)
        exp_output_9 = {'outputDir': cwd,
                        'batchRun': False,
                        'overwrite': False,
                        'outFiles': 'all',
                        'filter': False,
                        'temperature': 100,
                        'resolution': None,
                        'PDT': 7.0,
                        'windowSize': 0.02,
                        'HETATM': False,
                        'removeAtoms': [],
                        'addAtoms': [],
                        'highlightAtoms': [],
                        'createOrigpdb': False,
                        'createAUpdb': False,
                        'createUCpdb': False,
                        'createAUCpdb': False,
                        'createTApdb': False}
        self.assertDictEqual(act_output_9, exp_output_9)

        # Check temperature='none' is converted to temperature=None
        input_10 = ['temperature=none']
        act_output_10 = parse_input_file_arguments(input_10)
        exp_output_10 = {'outputDir': cwd,
                         'batchRun': False,
                         'overwrite': False,
                         'outFiles': 'all',
                         'filter': False,
                         'temperature': None,
                         'resolution': None,
                         'PDT': 7.0,
                         'windowSize': 0.02,
                         'HETATM': False,
                         'removeAtoms': [],
                         'addAtoms': [],
                         'highlightAtoms': [],
                         'createOrigpdb': False,
                         'createAUpdb': False,
                         'createUCpdb': False,
                         'createAUCpdb': False,
                         'createTApdb': False}
        self.assertDictEqual(act_output_10, exp_output_10)

        # Check raises ArgumentError if value specified for resolution isn't
        # recognised
        input_11 = ['resolution=xyz']
        self.assertRaises(ArgumentError, parse_input_file_arguments, input_11)

        # Check resolution='none' is converted to temperature=None
        input_12 = ['resolution=none']
        act_output_12 = parse_input_file_arguments(input_12)
        exp_output_12 = {'outputDir': cwd,
                         'batchRun': False,
                         'overwrite': False,
                         'outFiles': 'all',
                         'filter': False,
                         'temperature': None,
                         'resolution': None,
                         'PDT': 7.0,
                         'windowSize': 0.02,
                         'HETATM': False,
                         'removeAtoms': [],
                         'addAtoms': [],
                         'highlightAtoms': [],
                         'createOrigpdb': False,
                         'createAUpdb': False,
                         'createUCpdb': False,
                         'createAUCpdb': False,
                         'createTApdb': False}
        self.assertDictEqual(act_output_12, exp_output_12)

        # Check raises ArgumentError if value specified for PDT isn't recognised
        input_13 = ['pdt=xyz']
        self.assertRaises(ArgumentError, parse_input_file_arguments, input_13)

        # Check raises ArgumentError if value specified for windowsize isn't
        # recognised
        input_14 = ['windowsize=xyz']
        self.assertRaises(ArgumentError, parse_input_file_arguments, input_14)

        # Check raises ArgumentError if value specified for windowsize isn't
        # a float in the range 0 < windowsize < 1
        input_15 = ['windowsize=1']
        self.assertRaises(ArgumentError, parse_input_file_arguments, input_15)

        # Check raises ArgumentError if value is specified for
        # proteinornucleicacid
        input_16 = ['proteinornucleicacid=protein']
        self.assertRaises(ArgumentError, parse_input_file_arguments, input_16)

        # Check raises ArgumentError if value specified for HETATM isn't
        # recognised (expect "keep" or "remove")
        input_17 = ['hetatm=True']
        self.assertRaises(ArgumentError, parse_input_file_arguments, input_17)

        # Check raises ArgumentError if value specified for removeatoms isn't
        # recognised
        input_18 = ['removeAtoms=A-C']
        self.assertRaises(ArgumentError, parse_input_file_arguments, input_18)

        # Check raises ArgumentError if value specified for addatoms isn't
        # recognised
        input_19 = ['addAtoms=1-2-3']
        self.assertRaises(ArgumentError, parse_input_file_arguments, input_19)

        # Check raises ArgumentError if value specified for highlightatoms isn't
        # recognised
        input_20 = ['highlightAtoms=1-A']
        self.assertRaises(ArgumentError, parse_input_file_arguments, input_20)

        # Check raises ArgumentError if value specified for createorigpdb isn't
        # recognised
        input_21 = ['createOrigpdb=Flase']
        self.assertRaises(ArgumentError, parse_input_file_arguments, input_21)

        # Check raises ArgumentError if value specified for createaupdb isn't
        # recognised
        input_22 = ['createAUpdb=Yse']
        self.assertRaises(ArgumentError, parse_input_file_arguments, input_22)

        # Check raises ArgumentError if value specified for createucpdb isn't
        # recognised
        input_23 = ['createUCpdb=1-6']
        self.assertRaises(ArgumentError, parse_input_file_arguments, input_23)

        # Check raises ArgumentError if value specified for createaucpdb isn't
        # recognised
        input_24 = ['createAUCpdb=blah']
        self.assertRaises(ArgumentError, parse_input_file_arguments, input_24)

        # Check raises ArgumentError if value specified for createtapdb isn't
        # recognised
        input_25 = ['createTApdb=np.inf']
        self.assertRaises(ArgumentError, parse_input_file_arguments, input_25)

        # Check raises ArgumentError if unrecognised variable provided
        input_26 = ['blah=1']
        self.assertRaises(ArgumentError, parse_input_file_arguments, input_26)
