
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


# python -m unittest tests/test_bnet_calculation.py

import unittest

from rabdam.Subroutines.CalculateBDamage import run_rabdam

class TestClass(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(TestClass, self).__init__(*args, **kwargs)
        self.pdbs = [
            '2O2X', '4EZF', '4MWU', '4MOV', '3NBM', '1GW1', '4EWE', '3F1P',
            '3IV0', '4ZWV', '1T2I', '3LX3', '5P4N', '5MAA', '1E73', '1YKI',
            '4WA4', '3V2J', '3CUI', '4XLA', '4DUK', '3V38', '1VJF', '5IO2',
            '5CM7', '2EHU', '5JOW', '2C54', '4GZK', '2NUM', '5FYO', '4ODK',
            '6EV4', '5P5U', '3VHV', '4JCK', '5EKM', '3H4O', '5JIG', '2H5S',
            '4M5I', '1Y59', '4C45', '5F90', '4NI3', '4Z6N', '5M2G', '5ER6',
            '4R0X', '5LLG', '1FCX', '5M90', '3NJK', '5CWG', '2P7O', '5SZC',
            '2I0K', '4RDK', '5MA0', '4C1E', '5EJ3', '2WUG', '4MPY', '4OTZ',
            '4IOO', '4Z6O', '4ZOT', '5PHB', '3UJC', '4FR8', '5PH8', '5UPM',
            '3MWX', '4KDX', '3WH5', '4P04', '5Y90', '4H31', '5HJE', '4YKK',
            '1Q0F', '5JP6', '1X7Y', '4ZC8', '5EPE', '4ZS9', '5VNX', '5IHV',
            '5J90', '4K6W', '3PBC', '5CMB', '4PSC', '5UPN', '4XLZ', '4XGY',
            '5M4G', '3A54'
        ]
        # Cif files and PDB files sometimes report resolution to different
        # numbers of significant figures => Bnet_percentile scores are slightly
        # different, hence the need to check both
        self.exp_bnet_cif_dict = {
            '2O2X': [3.300580966, 93.06184012],
            '4EZF': [3.193514624, 80.73899371],
            '4MWU': [3.185476349, 96.62425344],
            '4MOV': [3.144130191, 90.54054054],
            '3NBM': [3.141821366, 81.15318417],
            '1GW1': [3.089945766, 94.78186484],
            '4EWE': [3.08241654, 91.83823529],
            '3F1P': [3.060628186, 77.41071429],
            '3IV0': [3.054440912, 84.77443609],
            '4ZWV': [3.0181382, 89.94469583],
            '1T2I': [3.004830448, 73.07692308],
            '3LX3': [2.964902529, 90.1932713],
            '5P4N': [2.916582486, 76.51858568],
            '5MAA': [2.91219352, 94.72764346],
            '1E73': [2.850203561, 87.43086978],
            '1YKI': [2.797739814, 93.13690262],
            '4WA4': [2.720540993, 96.76280165],
            '3V2J': [2.669599635, 91.44807761],
            '3CUI': [2.666548749, 84.46455505],
            '4XLA': [2.618363208, 82.46135553],
            '4DUK': [2.84941657, 89.28928929],
            '3V38': [2.502413452, 80.69381599],
            '1VJF': [2.494688495, 88.9943074],
            '5IO2': [2.467587911, 58.45493562],
            '5CM7': [2.44869046, 81.6750179],
            '2EHU': [2.440416614, 90.80758245],
            '5JOW': [2.439619791, 84.03776978],
            '2C54': [2.379224017, 76.6214178],
            '4GZK': [2.349526276, 87.03757646],
            '2NUM': [2.659568414, 84.26344897],
            '5FYO': [2.319618192, 74.56008044],
            '4ODK': [2.304354685, 64.0459364],
            '6EV4': [2.302433369, 47.27954972],
            '5P5U': [2.288966997, 66.38554217],
            '3VHV': [2.287083451, 63.43984962],
            '4JCK': [2.27150332, 44.74370113],
            '5EKM': [2.258574341, 56.23978202],
            '3H4O': [2.231817033, 72.09653092],
            '5JIG': [2.247664542, 42.84591195],
            '2H5S': [2.206850226, 49.42630185],
            '4M5I': [2.169405117, 38.2751938],
            '1Y59': [2.138787261, 41.7167382],
            '4C45': [2.131256276, 65.94594595],
            '5F90': [2.111640041, 79.16666667],
            '4NI3': [2.088735516, 54.6819788],
            '4Z6N': [2.094327694, 68.30543933],
            '5M2G': [2.06566475, 83.56271098],
            '5ER6': [2.05707889, 69.50608447],
            '4R0X': [2.006996308, 34.67811159],
            '5LLG': [1.981501196, 29.00390625],
            '1FCX': [1.978251127, 61.2960761],
            '5M90': [1.964999541, 67.62589928],
            '3NJK': [1.87824308, 57.56661639],
            '5CWG': [1.949818624, 31.58798283],
            '2P7O': [1.921138477, 55.97484277],
            '5SZC': [1.962633169, 31.92090395],
            '2I0K': [1.901555841, 64.92805755],
            '4RDK': [1.886900766, 72.29601518],
            '5MA0': [1.877853781, 82.57178527],
            '4C1E': [1.877575448, 42.66784452],
            '5EJ3': [1.875439995, 31.23181377],
            '2WUG': [1.873476725, 76.78525058],
            '4MPY': [1.842338963, 80.68833652],
            '4OTZ': [1.813684419, 35.59168925],
            '4IOO': [1.779618661, 22.21214869],
            '4Z6O': [1.800528596, 67.9316888],
            '4ZOT': [1.799163077, 37.72084806],
            '5PHB': [1.783879628, 30.19417476],
            '3UJC': [1.747894856, 20.36163522],
            '4FR8': [1.738876799, 81.43302181],
            '5PH8': [1.737369291, 34.09893993],
            '5UPM': [1.736663507, 66.22349982],
            '3MWX': [1.732318778, 43.24324324],
            '4KDX': [1.726256735, 29.04135338],
            '3WH5': [1.749202392, 55.9352518],
            '4P04': [1.714107945, 79.39964685],
            '5Y90': [1.685599882, 20.48192771],
            '4H31': [1.676976403, 62.95364714],
            '5HJE': [1.664006923, 29.85865724],
            '4YKK': [1.654971189, 29.86270023],
            '1Q0F': [1.646880018, 78.56697819],
            '5JP6': [1.629246723, 42.58421317],
            '1X7Y': [1.64715517, 50.4004004],
            '4ZC8': [1.60606196, 64.49864499],
            '5EPE': [1.598585924, 70.21223471],
            '4ZS9': [1.582398487, 23.4375],
            '5VNX': [1.543824945, 50.0427716],
            '5IHV': [1.542271159, 10.41275797],
            '5J90': [1.526469901, 21.79585572],
            '4K6W': [1.520316883, 34.79135244],
            '3PBC': [1.512656216, 21.33867277],
            '5CMB': [1.50829636, 17.08737864],
            '4PSC': [1.491796934, 8.601216334],
            '5UPN': [1.477252783, 55.85562192],
            '4XLZ': [1.473298738, 33.29124579],
            '4XGY': [1.460754281, 30.65371025],
            '5M4G': [1.400380357, 25.39356605],
            '3A54': [1.32172782, 20.56309703]
        }

        self.exp_bnet_pdb_dict = {
            '2O2X': [3.300580966, 93.06184012],
            '4EZF': [3.193514624, 80.73899371],
            '4MWU': [3.185476349, 96.62425344],
            '4MOV': [3.144130191, 90.54054054],
            '3NBM': [3.141821366, 81.15318417],
            '1GW1': [3.089945766, 94.78186484],
            '4EWE': [3.08241654, 91.83823529],
            '3F1P': [3.060628186, 77.41071429],
            '3IV0': [3.054440912, 84.77443609],
            '4ZWV': [3.0181382, 89.94469583],
            '1T2I': [3.004830448, 73.07692308],
            '3LX3': [2.964902529, 90.1932713],
            '5P4N': [2.916582486, 76.51858568],
            '5MAA': [2.91219352, 94.72764346],
            '1E73': [2.850203561, 87.43086978],
            '1YKI': [2.797739814, 93.13690262],
            '4WA4': [2.720540993, 96.76280165],
            '3V2J': [2.669599635, 91.44807761],
            '3CUI': [2.666548749, 84.46455505],
            '4XLA': [2.618363208, 82.46135553],
            '4DUK': [2.84941657, 89.28928929],
            '3V38': [2.502413452, 80.69381599],
            '1VJF': [2.494688495, 88.9943074],
            '5IO2': [2.467587911, 58.45493562],
            '5CM7': [2.44869046, 81.6750179],
            '2EHU': [2.440416614, 90.80758245],
            '5JOW': [2.439619791, 84.03776978],
            '2C54': [2.379224017, 76.6214178],
            '4GZK': [2.349526276, 87.03757646],
            '2NUM': [2.659568414, 84.26344897],
            '5FYO': [2.319618192, 74.56008044],
            '4ODK': [2.304354685, 64.0459364],
            '6EV4': [2.302433369, 47.27954972],
            '5P5U': [2.288966997, 66.90255852],
            '3VHV': [2.287083451, 63.43984962],
            '4JCK': [2.27150332, 44.74370113],
            '5EKM': [2.258574341, 56.23978202],
            '3H4O': [2.231817033, 72.09653092],
            '5JIG': [2.247664542, 42.84591195],
            '2H5S': [2.206850226, 49.42630185],
            '4M5I': [2.169405117, 38.2751938],
            '1Y59': [2.138787261, 41.7167382],
            '4C45': [2.131256276, 65.94594595],
            '5F90': [2.111640041, 79.16666667],
            '4NI3': [2.088735516, 54.6819788],
            '4Z6N': [2.094327694, 68.30543933],
            '5M2G': [2.06566475, 83.56271098],
            '5ER6': [2.05707889, 69.50608447],
            '4R0X': [2.006996308, 34.67811159],
            '5LLG': [1.981501196, 29.00390625],
            '1FCX': [1.978251127, 61.2960761],
            '5M90': [1.964999541, 67.62589928],
            '3NJK': [1.87824308, 57.56661639],
            '5CWG': [1.949818624, 31.58798283],
            '2P7O': [1.921138477, 55.97484277],
            '5SZC': [1.962633169, 31.36792453],
            '2I0K': [1.901555841, 64.92805755],
            '4RDK': [1.886900766, 72.29601518],
            '5MA0': [1.877853781, 82.57178527],
            '4C1E': [1.877575448, 42.66784452],
            '5EJ3': [1.875439995, 31.81818182],
            '2WUG': [1.873476725, 76.78525058],
            '4MPY': [1.842338963, 80.68833652],
            '4OTZ': [1.813684419, 35.59168925],
            '4IOO': [1.779618661, 22.21214869],
            '4Z6O': [1.800528596, 67.9316888],
            '4ZOT': [1.799163077, 37.72084806],
            '5PHB': [1.783879628, 30.19417476],
            '3UJC': [1.747894856, 20.36163522],
            '4FR8': [1.738876799, 81.43302181],
            '5PH8': [1.737369291, 34.09893993],
            '5UPM': [1.736663507, 66.22349982],
            '3MWX': [1.732318778, 43.24324324],
            '4KDX': [1.726256735, 29.04135338],
            '3WH5': [1.749202392, 55.9352518],
            '4P04': [1.714107945, 79.39964685],
            '5Y90': [1.685599882, 20.48192771],
            '4H31': [1.676976403, 62.95364714],
            '5HJE': [1.664006923, 29.85865724],
            '4YKK': [1.654971189, 29.86270023],
            '1Q0F': [1.646880018, 78.56697819],
            '5JP6': [1.629246723, 42.58421317],
            '1X7Y': [1.64715517, 50.4004004],
            '4ZC8': [1.60606196, 64.49864499],
            '5EPE': [1.598585924, 70.21223471],
            '4ZS9': [1.582398487, 23.4375],
            '5VNX': [1.543824945, 50.0427716],
            '5IHV': [1.542271159, 10.41275797],
            '5J90': [1.526469901, 22.19931271],
            '4K6W': [1.520316883, 34.79135244],
            '3PBC': [1.512656216, 21.33867277],
            '5CMB': [1.50829636, 17.08737864],
            '4PSC': [1.491796934, 8.601216334],
            '5UPN': [1.477252783, 55.85562192],
            '4XLZ': [1.473298738, 33.29124579],
            '4XGY': [1.460754281, 30.38827258],
            '5M4G': [1.400380357, 25.39356605],
            '3A54': [1.32172782, 20.56309703]
        }

    def test_bnet_values(self):
        """
        Checks that RABDAM calculates expected Bnet and Bnet-percentile values
        for a selection of PDB entries
        """

        import os
        import requests
        import shutil

        if not os.path.isdir('tests/temp_files/'):
            os.mkdir('tests/temp_files/')

        for code in self.pdbs:
            # Checks cif file
            exp_bnet_cif = self.exp_bnet_cif_dict[code][0]
            exp_bnet_p_cif = self.exp_bnet_cif_dict[code][1]

            cif_text = requests.get('https://files.rcsb.org/view/%s.cif' % code)
            with open('tests/temp_files/%s.cif' % code, 'w') as f:
                f.write(cif_text.text)
            rabdam_cif = run_rabdam(
                pathToInput='%s/tests/temp_files/%s.cif' % (os.getcwd(), code),
                outputDir='%s/tests/temp_files/' % os.getcwd(),
                batchRun=True,
                overwrite=True,
                outFiles="all",
                filterInput=False,
                temperature=None,
                resolution=None,
                PDT=7,
                windowSize=0.02,
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
            success = rabdam_cif.rabdam_dataframe(test=True)
            (
                act_prot_bnet, act_prot_bnet_p, act_na_bnet, act_na_bnet_p
            ) = rabdam_cif.rabdam_analysis()

            self.assertTrue(success)
            self.assertEqual(round(exp_bnet_cif, 6), round(act_prot_bnet, 6))
            self.assertEqual(round(exp_bnet_p_cif, 6), round(act_prot_bnet_p, 6))
            self.assertIsNone(act_na_bnet)
            self.assertIsNone(act_na_bnet_p)
            os.remove('tests/temp_files/%s.cif' % code)
            os.remove('tests/temp_files/Logfiles/Bnet_protein.pkl')

            # Checks PDB file
            exp_bnet_pdb = self.exp_bnet_pdb_dict[code][0]
            exp_bnet_p_pdb = self.exp_bnet_pdb_dict[code][1]

            pdb_text = requests.get('https://files.rcsb.org/view/%s.pdb' % code)
            with open('tests/temp_files/%s.pdb' % code, 'w') as f:
                f.write(pdb_text.text)
            rabdam_pdb = run_rabdam(
                pathToInput='%s/tests/temp_files/%s.pdb' % (os.getcwd(), code),
                outputDir='%s/tests/temp_files/' % os.getcwd(),
                batchRun=True,
                overwrite=True,
                outFiles="all",
                filterInput=False,
                temperature=None,
                resolution=None,
                PDT=7,
                windowSize=0.02,
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
            success = rabdam_pdb.rabdam_dataframe(test=True)
            (
                act_prot_bnet, act_prot_bnet_p, act_na_bnet, act_na_bnet_p
            ) = rabdam_pdb.rabdam_analysis()

            self.assertTrue(success)
            self.assertEqual(round(exp_bnet_pdb, 6), round(act_prot_bnet, 6))
            self.assertEqual(round(exp_bnet_p_pdb, 6), round(act_prot_bnet_p, 6))
            self.assertIsNone(act_na_bnet)
            self.assertIsNone(act_na_bnet_p)
            os.remove('tests/temp_files/%s.pdb' % code)
            os.remove('tests/temp_files/Logfiles/Bnet_protein.pkl')

        shutil.rmtree('tests/temp_files/')
