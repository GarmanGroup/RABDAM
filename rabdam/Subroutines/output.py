
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


class generate_output_files(object):
    def __init__(self, pdb_file_path, df):
        self.pdb_file_path = pdb_file_path
        self.pdb_code = pdb_file_path.split('/')[-1]
        self.df = df

    def make_csv(self, bdamatomList, window):
        # Returns a csv file containing a complete set of atom information
        # (including both that provided in the input PDB / mmCif file and also
        # the BDamage values calculated by RABDAM) for all atoms considered for
        # BDamage analysis. (This provides the user with a copy of the raw data
        # which they can manipulate as they wish.)

        newFile = open('%s_BDamage.csv' % self.pdb_file_path, 'w')

        # Defines column header abbreviations
        newFile.write('REC = RECORD NAME\n'
                      'ATMNUM = ATOM SERIAL NUMBER\n'
                      'ATMNAME = ATOM NAME\n'
                      'CONFORMER = ALTERNATE LOCATION INDICATOR\n'
                      'RESNAME = RESIDUE NAME\n'
                      'CHAIN = CHAIN IDENTIFIER\n'
                      'RESNUM = RESIDUE SEQUENCE NUMBER\n'
                      'INSCODE = INSERTION CODE'
                      'XPOS = ORTHOGONAL COORDINATES FOR X IN ANGSTROMS\n'
                      'YPOS = ORTHOGONAL COORDINATES FOR Y IN ANGSTROMS\n'
                      'ZPOS = ORTHOGONAL COORDINATES FOR Z IN ANGSTROMS\n'
                      'OCC = OCCUPANCY\n'
                      'BFAC = B FACTOR (TEMPERATURE FACTOR)\n'
                      'ELEMENT = ELEMENT SYMBOL\n'
                      'CHARGE = CHARGE ON ATOM\n'
                      'PD = PACKING DENSITY (ATOMIC CONTACT NUMBER)\n')
        newFile.write('AVRG_BF = AVERAGE B FACTOR FOR ATOMS IN A SIMILAR '
                      'PACKING DENSITY ENVIRONMENT (SLIDING WINDOW SIZE '
                      '= %s)\n' % window)
        newFile.write('BDAM = B DAMAGE VALUE\n'
                      '\n')
        # Writes column headers
        newFile.write('REC' + ','
                      'ATMNUM' + ','
                      'ATMNAME' + ','
                      'CONFORMER' + ','
                      'RESNAME' + ','
                      'CHAIN' + ','
                      'RESNUM' + ','
                      'INSCODE' + ','
                      'XPOS' + ','
                      'YPOS' + ','
                      'ZPOS' + ','
                      'OCC' + ','
                      'BFAC' + ','
                      'ELEMENT' + ','
                      'CHARGE' + ','
                      'PD' + ','
                      'AVRG_BF' + ','
                      'BDAM' + '\n')

        # Writes properties of each atom considered for BDamage analysis.
        for atm in bdamatomList:
            newFile.write(atm.lineID + ',')
            newFile.write(str(atm.atomNum) + ',')
            newFile.write(atm.atomType + ',')
            newFile.write(atm.conformer + ',')
            newFile.write(atm.resiType + ',')
            newFile.write(atm.chainID + ',')
            newFile.write(str(atm.resiNum) + ',')
            newFile.write(atm.insCode + ',')
            newFile.write(str(atm.xyzCoords[0][0]) + ',')
            newFile.write(str(atm.xyzCoords[1][0]) + ',')
            newFile.write(str(atm.xyzCoords[2][0]) + ',')
            newFile.write(str(atm.occupancy) + ',')
            newFile.write(str(atm.bFactor) + ',')
            newFile.write(atm.atomID + ',')
            newFile.write(str(atm.charge) + ',')
            newFile.write(str(atm.pd) + ',')
            newFile.write(str(atm.avrg_bf) + ',')
            newFile.write(str(atm.bd) + '\n')

        newFile.close()

    def generate_cif_lines(self, cif_column_widths):
        cif_lines = []

        for row in range(self.df.shape[0]):
            lineID = self.df.REC[row].ljust(cif_column_widths['REC']) + ' '
            atomNum = str(self.df.ATMNUM[row]).ljust(cif_column_widths['ATMNUM']) + ' '
            atomType = self.df.ATMNAME[row].ljust(cif_column_widths['ATMNAME']) + ' '
            conformer = self.df.CONFORMER[row].ljust(cif_column_widths['CONFORMER']) + ' '
            resiType = self.df.RESNAME[row].ljust(cif_column_widths['RESNAME']) + ' '
            chainID = self.df.CHAIN[row].ljust(cif_column_widths['CHAIN']) + ' '
            resiNum = str(self.df.RESNUM[row]).ljust(cif_column_widths['RESNUM']) + ' '
            insCode = self.df.INSCODE[row].ljust(cif_column_widths['INSCODE']) + ' '
            x_coord = str(self.df.XPOS[row]).ljust(cif_column_widths['XPOS']) + ' '
            y_coord = str(self.df.YPOS[row]).ljust(cif_column_widths['YPOS']) + ' '
            z_coord = str(self.df.ZPOS[row]).ljust(cif_column_widths['ZPOS']) + ' '
            occupancy = str(self.df.OCC[row]).ljust(cif_column_widths['OCC']) + ' '
            bFactor = str(self.df.BFAC[row]).ljust(cif_column_widths['BFAC']) + ' '
            bDamage = '{0:.2f}'.format(self.df.BDAM[row]).ljust(cif_column_widths['BDAM']) + ' '
            element = self.df.ELEMENT[row].ljust(cif_column_widths['ELEMENT']) + ' '
            charge = self.df.CHARGE[row].ljust(cif_column_widths['CHARGE'])

            cif_line = (lineID + atomNum + atomType + conformer + resiType
                        + chainID + resiNum + insCode + x_coord + y_coord
                        + z_coord + occupancy + bFactor + bDamage + element
                        + charge)
            cif_lines.append(cif_line)

        return cif_lines

    def generate_aniso_cif_lines(self, aniso_rec):
        skip = False

        if aniso_rec:
            aniso_lines = ['#', 'loop_']

            try:
                columns = [line.split('.')[1] for line in aniso_rec if
                           line.startswith('_atom_site_anisotrop')]
                atomNum_index = columns.index('id')
            except ValueError:
                skip = True
                print('\n\nERROR: mmCIF file _atom_site_anisotrop labels do '
                      'not follow the expected format\n. Skipping cif file '
                      'output\n.')

            for label in columns:
                label = '_atom_site_anisotrop.' + label
                aniso_lines.append(label)

            for line in aniso_rec:
                if not line[0:5] in ['loop_', '_atom']:
                    values = line.split()
                    if int(values[atomNum_index]) in self.df.ATMNUM.tolist():
                        aniso_lines.append(line)
        else:
            aniso_lines = []

        return aniso_lines, skip

    def write_output_cif(self, cif_header_lines, cif_lines, aniso_lines,
                         cif_footer_lines):
        # Writes an mmCif file for the input structure with an additional
        # column of BDamage values for those atoms included in the calculation.
        new_cif = open('%s_BDamage.cif' % self.pdb_file_path, 'w')

        for line in cif_header_lines:
            if line.replace(' ', '') != '':
                new_cif.write('%s\n' % line)

        new_cif.write('loop_\n')

        # Writes column labels to output cif file
        new_cif.write('_atom_site.group_PDB\n' +
                      '_atom_site.id\n' +
                      '_atom_site.auth_atom_id\n' +
                      '_atom_site.label_alt_id\n' +
                      '_atom_site.auth_comp_id\n' +
                      '_atom_site.auth_asym_id\n' +
                      '_atom_site.auth_seq_id\n' +
                      '_atom_site.pdbx_PDB_ins_code\n' +
                      '_atom_site.Cartn_x\n' +
                      '_atom_site.Cartn_y\n' +
                      '_atom_site.Cartn_z\n' +
                      '_atom_site.occupancy\n' +
                      '_atom_site.B_iso_or_equiv\n' +
                      '_atom_site.B_Damage\n' +
                      '_atom_site.type_symbol\n' +
                      '_atom_site.pdbx_formal_charge\n')

        # Writes ATOM / HETATM records to output cif file
        for line in cif_lines:
            new_cif.write('%s\n' % line)

        for line in aniso_lines:
            new_cif.write('%s\n' % line)

        for line in cif_footer_lines:
            new_cif.write('%s\n' % line)

        new_cif.close()

    def make_histogram(self, highlightAtoms):
        # Returns a kernel density estimate of the BDamage values of every atom
        # considered for BDamage analysis. Any atom whose number is listed
        # in the highlightAtoms option in the input file will be marked on the
        # plot. (Note that it is recommended no more than 6 atoms are listed
        # in the highlightAtoms option in the input file (beyond 6 atoms, the
        # colour scheme will repeat itself, and in addition the key may not fit
        # onto the graph).)

        import sys
        # Ensures that seaborn uses scipy for kde bandwidth calculation rather
        # than statsmodels
        sys.modules['statsmodels']=None
        import matplotlib.pyplot as plt
        import seaborn as sns

        # Sets figure aesthetics to seaborn defaults
        sns.set()

        plt.clf()  # Prevents the kernel density estimate of all atoms
        # considered for BDamage analysis from being plotted on the same axes
        # as the kernel density estimate of the atoms considered for
        # calculation of the Bnet summary metric

        # Generates kernel density plot
        plot = sns.distplot(self.df.BDAM.values, hist=False, rug=True,
                            kde_kws={'bw':'scott'})

        # Marks on the positions of any atoms whose numbers are listed in the
        # highlightAtoms option specified in the input file.
        xy_values = plot.get_lines()[0].get_data()
        y_values = xy_values[1]

        highlighted_atoms = []
        for number in highlightAtoms:
            for index, value in enumerate(self.df.ATMNUM.values):
                if int(number) == value:
                    b_dam_value = self.df.BDAM.values[index]
                    line, = plt.plot([b_dam_value, b_dam_value],
                                     [0, max(y_values)], linewidth=2,
                                     label=' atom ' + str(number) +
                                     '\n BDamage = {:.2f}'.format(b_dam_value))
                    highlighted_atoms.append(line)
                    break

        if len(highlighted_atoms) >= 1:
            plt.legend(handles=highlighted_atoms)
        plt.xlabel('B Damage')
        plt.ylabel('Normalised Frequency')
        plt.title(self.pdb_code + ' BDamage kernel density plot')
        plt.savefig(self.pdb_file_path + '_BDamage.svg')

    def calculate_Bnet(self, window_name, pdt_name, window):
        # Plots a kernel density estimate of the BDamage values of Glu O and
        # Asp O atoms. The summary metric Bnet is then calculated as the ratio
        # of the areas under the curve either side of the median (of the
        # overall BDamage distribution).

        import os
        import sys
        # Ensures that seaborn uses scipy for kde bandwidth calculation rather
        # than statsmodels
        sys.modules['statsmodels']=None
        import pandas as pd
        import matplotlib.pyplot as plt
        import seaborn as sns

        # Sets figure aesthetics to seaborn defaults
        sns.set()

        # Selects Glu / Asp terminal oxygen atoms from complete DataFrame.
        prot = self.df[(self.df.RESNAME.isin(['GLU', 'ASP']))
                       & (self.df.ATMNAME.isin(['OE1', 'OE2', 'OD1', 'OD2']))]
        # Selects atoms of sugar-phosphate C-O bonds from complete DataFrame.
        na = self.df[self.df.ATMNAME.isin(["O3'", "O5'", "C3'", "C5'"])]

        if prot.empty and na.empty:
            print('\n\nERROR: No sites used for Bnet calculation present in structure\n')
            return

        if not prot.empty:
            # Calculates median of protein BDamage distribution
            median = self.df.BDAM.median()

            plt.clf()  # Prevents the kernel density estimate of the atoms
            # considered for calculation of the Bnet summary metric from being
            # plotted on the same axes as the kernel density estimate of all
            # atoms considered for BDamage analysis.
            plot = sns.distplot(prot.BDAM.values, hist=False, rug=True,
                                kde_kws={'bw':'scott'})
            plt.xlabel('B Damage')
            plt.ylabel('Normalised Frequency')
            plt.title(self.pdb_code + ' Bnet kernel density plot')

            # Extracts an array of 128 (x, y) coordinate pairs evenly spaced
            # along the x(BDamage)-axis from the kernel density plot. These
            # coordinate pairs are used to calculate, via the trapezium rule,
            # the area under the curve between the smallest value of x and the
            # median (= total_area_LHS), and the area under the curve between
            # the median and the largest value of x (= total_area_RHS). The
            # Bnet summary metric is then calculated as the ratio of
            # total_area_RHS to total_area_LHS.
            xy_values = plot.get_lines()[0].get_data()
            x_values = xy_values[0]
            y_values = xy_values[1]
            height = (x_values[-1]-x_values[0]) / (len(x_values)-1)

            total_area_LHS = 0
            for index, value in enumerate(y_values):
                if x_values[index] < median:
                    area_LHS = (((y_values[index] + y_values[index+1]) / 2)
                                * height)
                    total_area_LHS = total_area_LHS + area_LHS

            total_area_RHS = 0
            for index, value in enumerate(y_values):
                if x_values[index] >= median and index < len(y_values)-1:
                    area_RHS = (((y_values[index] + y_values[index+1]) / 2)
                                * height)
                    total_area_RHS = total_area_RHS + area_RHS

            # Calculates area ratio (= Bnet)
            ratio = total_area_RHS / total_area_LHS

            plt.annotate('Bnet = {:.1f}'.format(ratio),
                         xy=(max(x_values)*0.65, max(y_values)*0.9),
                         fontsize=10)
            plt.annotate('Median = {:.2f}'.format(median),
                         xy=(max(x_values)*0.65, max(y_values)*0.85),
                         fontsize=10)
            plt.savefig(self.pdb_file_path + '_Bnet_Protein.svg')

            if not os.path.isfile('Logfiles/Bnet_Protein.csv'):
                Bnet_list = open('Logfiles/Bnet_Protein.csv', 'w')
                Bnet_list.write('PDB' + ',')
                Bnet_list.write('Bnet' + ',')
                Bnet_list.write('Window_size' + ',')
                Bnet_list.write('Window_size (%)' + ',')
                Bnet_list.write('PDT' + '\n')
                Bnet_list.close()
            Bnet_list = open('Logfiles/Bnet_Protein.csv', 'a')
            Bnet_list.write('%s' % self.pdb_code + ',')
            Bnet_list.write('%s' % ratio + ',')
            Bnet_list.write('%s' % window + ',')
            Bnet_list.write('%s' % window_name + ',')
            Bnet_list.write('%s' % pdt_name + '\n')
            Bnet_list.close()

        if not na.empty:
            # Calculates median of nucleic acid BDamage distribution
            median = self.na.BDAM.median()

            plt.clf()  # Prevents the kernel density estimate of the atoms
            # considered for calculation of the Bnet summary metric from being
            # plotted on the same axes as the kernel density estimate of all
            # atoms considered for BDamage analysis.
            plot = sns.distplot(na.BDAM.values, hist=False, rug=True,
                                kde_kws={'bw':'scott'})
            plt.xlabel('B Damage')
            plt.ylabel('Normalised Frequency')
            plt.title(self.pdb_code + ' Bnet kernel density plot')

            # Extracts an array of 128 (x, y) coordinate pairs evenly spaced
            # along the x(BDamage)-axis from the kernel density plot. These
            # coordinate pairs are used to calculate, via the trapezium rule,
            # the area under the curve between the smallest value of x and the
            # median (= total_area_LHS), and the area under the curve between
            # the median and the largest value of x (= total_area_RHS). The
            # Bnet summary metric is then calculated as the ratio of
            # total_area_RHS to total_area_LHS.
            xy_values = plot.get_lines()[0].get_data()
            x_values = xy_values[0]
            y_values = xy_values[1]
            height = (x_values[-1]-x_values[0]) / (len(x_values)-1)

            total_area_LHS = 0
            for index, value in enumerate(y_values):
                if x_values[index] < median:
                    area_LHS = (((y_values[index] + y_values[index+1]) / 2)
                                * height)
                    total_area_LHS = total_area_LHS + area_LHS

            total_area_RHS = 0
            for index, value in enumerate(y_values):
                if x_values[index] >= median and index < len(y_values)-1:
                    area_RHS = (((y_values[index] + y_values[index+1]) / 2)
                                * height)
                    total_area_RHS = total_area_RHS + area_RHS

            # Calculates area ratio (= Bnet)
            ratio = total_area_RHS / total_area_LHS

            plt.annotate('Bnet = {:.1f}'.format(ratio),
                         xy=(max(x_values)*0.65, max(y_values)*0.9),
                         fontsize=10)
            plt.annotate('Median = {:.2f}'.format(median),
                         xy=(max(x_values)*0.65, max(y_values)*0.85),
                         fontsize=10)
            plt.savefig(str(self.pdb_file_path)+'_Bnet_NA.svg')

            if not os.path.isfile('Logfiles/Bnet_NA.csv'):
                Bnet_list = open('Logfiles/Bnet_NA.csv', 'w')
                Bnet_list.write('PDB' + ',')
                Bnet_list.write('Bnet' + ',')
                Bnet_list.write('Window_size' + ',')
                Bnet_list.write('Window_size (%)' + ',')
                Bnet_list.write('PDT' + '\n')
                Bnet_list.close()
            Bnet_list = open('Logfiles/Bnet_NA.csv', 'a')
            Bnet_list.write('%s' % self.pdb_code + ',')
            Bnet_list.write('%s' % ratio + ',')
            Bnet_list.write('%s' % window + ',')
            Bnet_list.write('%s' % window_name + ',')
            Bnet_list.write('%s' % pdt_name + '\n')
            Bnet_list.close()

    def write_html_summary(self, output_options, highlightAtoms):
        # Writes an html file summarising the results generated by RABDAM.
        import os
        import time
        import pandas as pd
        import numpy as np

        # Records time (for reporting in summary html file)
        tm = time.localtime()
        day = tm.tm_mday
        month = tm.tm_mon
        year = tm.tm_year
        hour = tm.tm_hour
        mins = tm.tm_min
        secs = tm.tm_sec

        # Filters the complete DataFrame to retain only the atoms specified to
        # be highlighted (default = none) by the user in the input file
        if len(highlightAtoms) != 0:
            highlightAtoms = [int(number) for number in highlightAtoms]
            sub_df_highlight = self.df[self.df.ATMNUM.isin(highlightAtoms)]
            sub_df_highlight = sub_df_highlight.round({'AVRG BF': 2, 'BDAM': 2})

        # Filter the complete DataFrame to retain the atoms with highest
        # BDamage values (> 1.96o = 2-tailed 95% confidence interval)
        sorted_df = self.df.sort_values(by='BDAM', ascending=False)
        b_dam_values = sorted_df['BDAM']
        ln_b_dam_values = np.log(b_dam_values)
        mean = np.mean(ln_b_dam_values)
        std_dev = np.std(ln_b_dam_values)
        cut = mean + (1.96*std_dev)
        cut_index = 0
        for index, value in enumerate(ln_b_dam_values.tolist()):
            if value > cut:
                cut_index = index
            else:
                break
        sub_df_top_site = sorted_df.drop(sorted_df.index[cut_index+1:])
        sub_df_top_site = sub_df_top_site.round({'AVRG BF': 2, 'BDAM': 2})

        # Filters the complete DataFrame to retain only Glu and Asp terminal
        # oxygen atoms.
        sub_df_prot = self.df[(self.df.RESNAME.isin(['GLU', 'ASP']))
                              & (self.df.ATMNAME.isin(['OE1', 'OE2', 'OD1', 'OD2']))]
        sub_df_prot = sub_df_prot.round({'AVRG BF': 2, 'BDAM': 2})
        # Filters the complete DataFrame to retain only atoms of
        # sugar-phosphate C-O bonds
        sub_df_na = self.df[self.df.ATMNAME.isin(["O3'", "O5'", "C3'", "C5'"])]
        sub_df_na = sub_df_na.round({'AVRG BF': 2, 'BDAM': 2})

        # Opens summary html file
        html_file = open(self.pdb_file_path + '_BDamage.html', 'w')
        # Specifies html file header information
        html_file.write('<!DOCTYPE html>\n'
                        '<html>\n'
                        '  <head>\n'
                        '    <script type="text/javascript" src="http://ajax.googleapis.com/ajax/libs/jquery/1.10.2/jquery.min.js"></script>\n'
                        '    <script type="text/javascript" src="file://%s/HTML_stylesheet.js"></script>\n' % os.path.dirname(os.path.abspath(__file__)))  # Locates
                                                                                                              # HTML stylesheets when RABDAM is run either as a
                                                                                                              # package or as a script
        html_file.write('    <link href="file://%s/HTML_stylesheet.css" type="text/css" rel="stylesheet">\n' % os.path.dirname(os.path.abspath(__file__)))
        html_file.write('    <title>'+self.pdb_code.replace('_', ' ')+' BDamage summary file</title>\n'
                        '  </head>\n')
        # Writes html file title
        html_file.write('  <body>\n'
                        '    <div class="file_name">\n'
                        '      <h1>'+self.pdb_code.replace('_', ' ')+' B<sub>Damage</sub> summary file</h1>\n')
        html_file.write('      <p id="file_info">Created on %02.0f/%02.0f/%04.0f at %02.0f:%02.0f:%02.0f.</p>\n'
                        % (day, month, year, hour, mins, secs))
        if 'csv' in output_options:
            html_file.write('      <p id="file_info">**NOTE: Table header abbreviations are as defined in '
                            '<a href="%s_BDamage.csv">%s_Bdamage.csv</a>.**</p>\n' % (self.pdb_code, self.pdb_code))
        html_file.write('    </div>\n')
        # Writes html file BDamage summary
        html_file.write('    <div>\n'
                        '      <ul class="main_list">\n'
                        '        <li><h2>B<sub>Damage</sub> distribution</h2></li>\n'
                        '        <li><img class="graphs" src="%s_BDamage.svg" /></li>\n' % self.pdb_code)
        if len(highlightAtoms) != 0:
            html_file.write('        <div>\n'
                            '          <ul class="sublist">\n'
                            '            <li><h3>Properties of highlighted atoms</h3></li>\n'
                            '            <li>'+sub_df_highlight.to_html(index=False)+'</li></br>\n'
                            '          </ul>\n'
                            '        </div>\n')
        html_file.write('        <div>\n'
                        '          <ul class="sublist">\n'
                        '            <li><h3>Properties of sites with B<sub>Damage</sub> values larger than the 95% confidence limit (2-tailed) threshold</h3></li>\n'
                        '            <li>'+sub_df_top_site.to_html(index=False)+'</li>\n'
                        '          </ul>\n'
                        '        </div>\n'
                        '      </ul>\n'
                        '    </div>\n')
        # Writes html file Bnet (protein) summary
        if os.path.isfile('%s_Bnet_Protein.svg' % self.pdb_file_path):
            html_file.write('    <div>\n'
                            '      <ul class="main_list">\n'
                            '        <li><h2>B<sub>net</sub> distribution (protein)</h2></li>\n'
                            '        <li><img class="graphs" src="%s_Bnet_Protein.svg" /></li>\n' % self.pdb_code)
            html_file.write('        <ul class="sublist">\n'
                            '          <li><h3>Properties of atoms used in calculation of B<sub>net</sub> (protein) metric</h3></li>\n'
                            '          <li>'+sub_df_prot.to_html(index=False)+'</li>\n'
                            '        </ul>\n'
                            '      </ul>\n'
                            '    </div>\n')
        # Writes html file Bnet (nucleic acid) summary
        if os.path.isfile('%s_Bnet_NA.svg' % self.pdb_file_path):
            html_file.write('    <div>\n'
                            '      <ul class="main_list">\n'
                            '        <li><h2>B<sub>net</sub> distribution (nucleic acid)</h2></li>\n'
                            '        <li><img class="graphs" src="%s_Bnet_NA.svg" /></li>\n' % self.pdb_code)
            html_file.write('        <ul class="sublist">\n'
                            '          <li><h3>Properties of atoms used in calculation of B<sub>net</sub> (nucleic acid) metric</h3></li>\n'
                            '          <li>'+sub_df_na.to_html(index=False)+'</li>\n'
                            '        </ul>\n'
                            '      </ul>\n'
                            '    </div>\n')
        # Concludes and closes summary html file
        html_file.write('  </body>\n'
                        '</html>\n')
        html_file.close()
