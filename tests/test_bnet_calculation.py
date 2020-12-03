
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

# python -m unittest tests/test_bnet_calculation.py

import os
import unittest

from rabdam.Subroutines.CalculateBDamage import rabdam

class TestClass(unittest.TestCase):

    def test_bnet_values(self):
        """
        Checks that RABDAM calculates expected Bnet values for a selection of
        PDB entries
        """

        import os
        import requests
        import shutil
        import pandas as pd

        exp_bnet_dict = {'2O2X': 3.300580966,
                         '4EZF': 3.193514624,
                         '4MWU': 3.185476349,
                         '4MOV': 3.144130191,
                         '3NBM': 3.141821366,
                         '1GW1': 3.105626889,
                         '4EWE': 3.08241654,
                         '3F1P': 3.060628186,
                         '3IV0': 3.054440912,
                         '4ZWV': 3.017330004,
                         '1T2I': 3.004830448,
                         '3LX3': 2.962424378,
                         '5P4N': 2.916582486,
                         '5MAA': 2.91219352,
                         '1E73': 2.850203561,
                         '1YKI': 2.797739814,
                         '4WA4': 2.720540993,
                         '3V2J': 2.669599635,
                         '3CUI': 2.666605946,
                         '4XLA': 2.624366813,
                         '4DUK': 2.854175949,
                         '3V38': 2.500984382,
                         '1VJF': 2.496374854,
                         '5IO2': 2.467587911,
                         '5CM7': 2.44869046,
                         '2EHU': 2.448290431,
                         '5JOW': 2.439619791,
                         '2C54': 2.379224017,
                         '4GZK': 2.349526276,
                         '2NUM': 2.326904729,
                         '5FYO': 2.319618192,
                         '4ODK': 2.304354685,
                         '6EV4': 2.302433369,
                         '5P5U': 2.288966997,
                         '3VHV': 2.285877338,
                         '4JCK': 2.27150332,
                         '5EKM': 2.258574341,
                         '3H4O': 2.231817033,
                         '5JIG': 2.247664542,
                         '2H5S': 2.206850226,
                         '4M5I': 2.169405117,
                         '1Y59': 2.138787261,
                         '4C45': 2.131256276,
                         '5F90': 2.11287042,
                         '4NI3': 2.088735516,
                         '4Z6N': 2.083743584,
                         '5M2G': 2.06566475,
                         '5ER6': 2.05707889,
                         '4R0X': 2.006996308,
                         '5LLG': 1.981501196,
                         '1FCX': 1.976990791,
                         '5M90': 1.96542442,
                         '3NJK': 1.955577757,
                         '5CWG': 1.949818624,
                         '2P7O': 1.921138477,
                         '5SZC': 1.962633169,
                         '2I0K': 1.901555841,
                         '4RDK': 1.886900766,
                         '5MA0': 1.877853781,
                         '4C1E': 1.877575448,
                         '5EJ3': 1.875439995,
                         '2WUG': 1.87334953,
                         '4MPY': 1.842338963,
                         '4OTZ': 1.835716553,
                         '4IOO': 1.828349113,
                         '4Z6O': 1.800528596,
                         '4ZOT': 1.799163077,
                         '5PHB': 1.783879628,
                         '3UJC': 1.747894856,
                         '4FR8': 1.738876799,
                         '5PH8': 1.736825591,
                         '5UPM': 1.736663507,
                         '3MWX': 1.733132746,
                         '4KDX': 1.729650659,
                         '3WH5': 1.717975404,
                         '4P04': 1.714107945,
                         '5Y90': 1.695283923,
                         '4H31': 1.674014779,
                         '5HJE': 1.662869176,
                         '4YKK': 1.653894709,
                         '1Q0F': 1.646880018,
                         '5JP6': 1.629246723,
                         '1X7Y': 1.618817315,
                         '4ZC8': 1.60606196,
                         '5EPE': 1.604407869,
                         '4ZS9': 1.582398487,
                         '5VNX': 1.543824945,
                         '5IHV': 1.542271159,
                         '5J90': 1.526469901,
                         '4K6W': 1.520316883,
                         '3PBC': 1.512738972,
                         '5CMB': 1.504620762,
                         '4PSC': 1.491796934,
                         '5UPN': 1.477252783,
                         '4XLZ': 1.473298738,
                         '4XGY': 1.465885549,
                         '5M4G': 1.400219288,
                         '3A54': 1.319587779}

        if not os.path.isdir('tests/temp_files/'):
            os.mkdir('tests/temp_files/')

        for code, exp_bnet in exp_bnet_dict.items():
            # Checks cif file
            cif_text = requests.get('https://files.rcsb.org/view/%s.cif' % code)
            with open('tests/temp_files/%s.cif' % code, 'w') as f:
                f.write(cif_text.text)
            rabdam_run = rabdam(
                pathToInput='%s/tests/temp_files/%s.cif' % (os.getcwd(), code),
                outputDir='%s/tests/temp_files/' % os.getcwd(),
                batchRun=True,
                overwrite=True,
                PDT=7,
                windowSize=0.02,
                protOrNA='protein',
                HETATM=False,
                removeAtoms=[],
                addAtoms=[],
                highlightAtoms=[],
                createOrigpdb=False,
                createAUpdb=False,
                createUCpdb=False,
                createAUCpdb=False,
                createTApdb=False
            )
            rabdam_run.rabdam_dataframe(test=True)
            rabdam_run.rabdam_analysis(
                output_options=['csv', 'pdb', 'cif', 'kde', 'bnet', 'summary']
            )
            bnet_df = pd.read_pickle('tests/temp_files/Logfiles/Bnet_protein.pkl')
            act_bnet_cif = bnet_df['Bnet'].tolist()[-1]
            self.assertEqual(round(exp_bnet, 7), round(act_bnet_cif, 7))
            os.remove('tests/temp_files/%s.cif' % code)
            os.remove('tests/temp_files/Logfiles/Bnet_protein.pkl')

            # Checks PDB file
            pdb_text = requests.get('https://files.rcsb.org/view/%s.pdb' % code)
            with open('tests/temp_files/%s.pdb' % code, 'w') as f:
                f.write(pdb_text.text)
            rabdam_run = rabdam(
                pathToInput='%s/tests/temp_files/%s.pdb' % (os.getcwd(), code),
                outputDir='%s/tests/temp_files/' % os.getcwd(),
                batchRun=True,
                overwrite=True,
                PDT=7,
                windowSize=0.02,
                protOrNA='protein',
                HETATM=False,
                removeAtoms=[],
                addAtoms=[],
                highlightAtoms=[],
                createOrigpdb=False,
                createAUpdb=False,
                createUCpdb=False,
                createAUCpdb=False,
                createTApdb=False
            )
            rabdam_run.rabdam_dataframe(test=True)
            rabdam_run.rabdam_analysis(
                output_options=['csv', 'pdb', 'cif', 'kde', 'bnet', 'summary']
            )
            bnet_df = pd.read_pickle(
                '%s/tests/temp_files/Logfiles/Bnet_protein.pkl' % os.getcwd()
            )
            act_bnet_pdb = bnet_df['Bnet'].tolist()[-1]
            self.assertEqual(round(exp_bnet, 7), round(act_bnet_pdb, 7))
            os.remove('tests/temp_files/%s.pdb' % code)
            os.remove('tests/temp_files/Logfiles/Bnet_protein.pkl')

        shutil.rmtree('tests/temp_files/')
