
# RABDAM
# Copyright (C) 2017 Garman Group, University of Oxford

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


class generate_output_files():
    def __init__(self, pdb_file_path, df):
        self.pdb_file_path = pdb_file_path
        self.pdb_code = pdb_file_path.split('/')[-1]
        self.df = df

    def make_csv(self, bdamatomList, window):
        # Returns a csv file containing a complete set of atom information
        # (including both that provided in the input PDB file and also the
        # BDamage values calculated by RABDAM) for all atoms considered for
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

    def make_histogram(self, highlightAtoms):
        # Returns a kernel density estimate of the BDamage values of every atom
        # considered for BDamage analysis. Any atom whose number is listed
        # in the highlightAtoms option in the input file will be marked on the
        # plot. (Note that it is recommended no more than 6 atoms are listed
        # in the highlightAtoms option in the input file (beyond 6 atoms, the
        # colour scheme will repeat itself, and in addition the key may not fit
        # onto the graph).)

        import matplotlib.pyplot as plt
        import seaborn as sns
        # Sets figure aesthetics to seaborn defaults
        sns.set()

        plt.clf()  # Prevents the kernel density estimate of all atoms
        # considered for BDamage analysis from being plotted on the same axes
        # as the kernel density estimate of the atoms considered for
        # calculation of the Bnet summary metric

        # Generates kernel density plot
        plot = sns.distplot(self.df.BDAM.values, hist=False, rug=True)

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

    def calculate_Bnet(self, window_name, pdt_name, count, window):
        # Plots a kernel density estimate of the BDamage values of Glu O and
        # Asp O atoms. The summary metric Bnet is then calculated as the ratio
        # of the areas under the curve either side of the median (of the
        # overall BDamage distribution).

        import os
        import pandas as pd
        import matplotlib.pyplot as plt
        import seaborn as sns
        # Sets figure aesthetics to seaborn defaults
        sns.set()

        # Selects Glu / Asp terminal oxygen atoms from complete DataFrame.
        a = self.df[(self.df.RESNAME.isin(['GLU']))
                    & (self.df.ATMNAME.isin(['OE1', 'OE2']))]
        b = self.df[(self.df.RESNAME.isin(['ASP']))
                    & (self.df.ATMNAME.isin(['OD1', 'OD2']))]
        dataframes = [a, b]
        prot = pd.concat(dataframes)
        # Selects atoms of sugar-phosphate C-O bonds from complete DataFrame.
        na = self.df[self.df.ATMNAME.isin(["O3'", "O5'", "C3'", "C5'"])]

        if prot.empty and na.empty:
            print('\nNo sites used for Bnet calculation present in structure\n')

        if not prot.empty:
            # Calculates median of protein BDamage distribution
            median = self.df.BDAM.median()

            plt.clf()  # Prevents the kernel density estimate of the atoms
            # considered for calculation of the Bnet summary metric from being
            # plotted on the same axes as the kernel density estimate of all
            # atoms considered for BDamage analysis.
            plot = sns.distplot(prot.BDAM.values, hist=False, rug=True)
            plt.xlabel('B Damage')
            plt.ylabel('Normalised Frequency')
            plt.title(self.pdb_code + ' Bnet kernel density plot')

            # Extracts an array of 128 (x, y) coordinate pairs evenly spaced
            # along the x(BDamage)-axis from the kernel density plot. These
            # coordinate pairs are used to calculate, via the trapezium rule,
            # the area under the curve between the smallest value of x and the
            # median (= area LHS), and the area under the curve between the
            # median and the largest value of x (= area RHS). The Bnet summary
            # metric is then calculated as the ratio of area RHS to area LHS.
            xy_values = plot.get_lines()[0].get_data()
            x_values = xy_values[0]
            y_values = xy_values[1]

            total_area_LHS = 0
            for index, value in enumerate(y_values):
                if x_values[index] < median:
                    area_LHS = (((y_values[index] + y_values[index+1]) / 2)
                                * ((x_values[-1]-x_values[0]) / (len(x_values)-1)))
                    total_area_LHS = total_area_LHS + area_LHS

            total_area_RHS = 0
            for index, value in enumerate(y_values):
                if x_values[index] >= median and index < len(y_values)-1:
                    area_RHS = (((y_values[index] + y_values[index+1]) / 2)
                                * ((x_values[-1]-x_values[0]) / (len(x_values)-1)))
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

            if count > 1:
                if not os.path.isfile('Logfiles/Bnet_Protein.csv'):
                    Bnet_list = open('Logfiles/Bnet_Protein.csv', 'w')
                    Bnet_list.write('PDB' + ',')
                    Bnet_list.write('Bnet' + ',')
                    Bnet_list.write('Window_size' + ',')
                    Bnet_list.write('Window_size (%)' + ',')
                    Bnet_list.write('PDT' + ',')
                    Bnet_list.close()
                Bnet_list = open('Logfiles/Bnet_Protein.csv', 'a')
                Bnet_list.write('\n%s' % self.pdb_code + ',')
                Bnet_list.write('%s' % ratio + ',')
                Bnet_list.write('%s' % window + ',')
                Bnet_list.write('%s' % window_name + ',')
                Bnet_list.write('%s' % pdt_name + ',')
                Bnet_list.close()

        if not na.empty:
            # Calculates median of nucleic acid BDamage distribution
            median = self.na.BDAM.median()

            plt.clf()  # Prevents the kernel density estimate of the atoms
            # considered for calculation of the Bnet summary metric from being
            # plotted on the same axes as the kernel density estimate of all
            # atoms considered for BDamage analysis.
            plot = sns.distplot(na.BDAM.values, hist=False, rug=True)
            plt.xlabel('B Damage')
            plt.ylabel('Normalised Frequency')
            plt.title(self.pdb_code + ' Bnet kernel density plot')

            # Extracts an array of 128 (x, y) coordinate pairs evenly spaced
            # along the x(BDamage)-axis from the kernel density plot. These
            # coordinate pairs are used to calculate, via the trapezium rule,
            # the area under the curve between the smallest value of x and the
            # median (= area LHS), and the area under the curve between the
            # median and the largest value of x (= area RHS). The Bnet summary
            # metric is then calculated as the ratio of area RHS to area LHS.
            xy_values = plot.get_lines()[0].get_data()
            x_values = xy_values[0]
            y_values = xy_values[1]

            total_area_LHS = 0
            for index, value in enumerate(y_values):
                if x_values[index] < median:
                    area_LHS = (((y_values[index] + y_values[index+1]) / 2)
                                * ((x_values[-1]-x_values[0]) / (len(x_values)-1)))
                    total_area_LHS = total_area_LHS + area_LHS

            total_area_RHS = 0
            for index, value in enumerate(y_values):
                if x_values[index] >= median and index < len(y_values)-1:
                    area_RHS = (((y_values[index] + y_values[index+1]) / 2)
                                * ((x_values[-1]-x_values[0]) / (len(x_values)-1)))
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

            if count > 1:
                if not os.path.isfile('Logfiles/Bnet_NA.csv'):
                    Bnet_list = open('Logfiles/Bnet_NA.csv', 'w')
                    Bnet_list.write('PDB' + ',')
                    Bnet_list.write('Bnet' + ',')
                    Bnet_list.write('Window_size' + ',')
                    Bnet_list.write('Window_size (%)' + ',')
                    Bnet_list.write('PDT' + ',')
                    Bnet_list.close()
                Bnet_list = open('Logfiles/Bnet_NA.csv', 'a')
                Bnet_list.write('\n%s' % self.pdb_code + ',')
                Bnet_list.write('%s' % ratio + ',')
                Bnet_list.write('%s' % window + ',')
                Bnet_list.write('%s' % window_name + ',')
                Bnet_list.write('%s' % pdt_name + ',')
                Bnet_list.close()

    def write_html_summary(self, cwd, output_options, highlightAtoms):
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
            sub_df_highlight = self.df[(self.df.ATMNUM.isin(highlightAtoms))]

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
        sub_df_top_site = sorted_df.drop(sorted_df.index[cut_index+1:])

        # Filters the complete DataFrame to retain only Glu and Asp terminal
        # oxygen atoms.
        a = self.df[(self.df.RESNAME.isin(['GLU']))
                    & (self.df.ATMNAME.isin(['OE1', 'OE2']))]
        b = self.df[(self.df.RESNAME.isin(['ASP']))
                    & (self.df.ATMNAME.isin(['OD1', 'OD2']))]
        dataframes = [a, b]
        sub_df_prot = pd.concat(dataframes)
        # Filters the complete DataFrame to retain only atoms of
        # sugar-phosphate C-O bonds
        sub_df_na = self.df[self.df.ATMNAME.isin(["O3'", "O5'", "C3'", "C5'"])]

        # Opens summary html file
        html_file = open(self.pdb_file_path + '_BDamage.html', 'w')
        # Specifies html file header information
        html_file.write('<!DOCTYPE html>\n'
                        '<html>\n'
                        '  <head>\n'
                        '    <script type="text/javascript" src="http://ajax.googleapis.com/ajax/libs/jquery/1.10.2/jquery.min.js"></script>\n'
                        '    <script type="text/javascript" src="file://%s/Subroutines/HTML_stylesheet.js"></script>\n' % cwd)
        html_file.write('    <link href="file://%s/Subroutines/HTML_stylesheet.css" type="text/css" rel="stylesheet">\n' % cwd)
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
