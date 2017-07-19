
class generate_output_files():
    def __init__(self, pdb_file_path, df):
        self.pdb_file_path = pdb_file_path
        self.pdb_code = pdb_file_path.split('/')[-1]
        self.df = df

    def make_csv(self, bdamatomList, window):
        # Returns a csv file containing a complete set of atom information
        # (including both that provided in the input PDB file and also the
        # B_damage values calculated by RABDAM) for all atoms considered for
        # B_damage analysis. (This provides the user with a copy of the raw data
        # which they can manipulate as they wish.)

        newFile = open('%sBdamage.csv' % self.pdb_file_path, 'w')

        # Defines column header abbreviations at top of file.
        newFile.write('REC = RECORD NAME\n'
                      'ATMNUM = ATOM SERIAL NUMBER\n'
                      'ATMNAME = ATOM NAME\n'
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
        newFile.write('AVRG_BF = AVERAGE B FACTOR FOR ATOMS IN A SIMILAR PACKING'
                      'DENSITY ENVIRONMENT (SLIDING WINDOW SIZE = %s)\n' % window)
        newFile.write('BDAM = BDAMAGE VALUE\n'
                      '\n')
        # Writes column headers to file
        newFile.write('REC' + ','
                      'ATMNUM' + ','
                      'ATMNAME' + ','
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

        # Writes properties of each atom considered fpr B_damage analysis to file.
        for atm in bdamatomList:
            newFile.write(str(atm.lineID) + ',')
            newFile.write(str(atm.atomNum) + ',')
            newFile.write(str(atm.atomType) + ',')
            newFile.write(str(atm.resiType) + ',')
            newFile.write(str(atm.chainID) + ',')
            newFile.write(str(atm.resiNum) + ',')
            newFile.write(str(atm.xyzCoords[0][0]) + ',')
            newFile.write(str(atm.xyzCoords[1][0]) + ',')
            newFile.write(str(atm.xyzCoords[2][0]) + ',')
            newFile.write(str(atm.occupancy) + ',')
            newFile.write(str(atm.bFactor) + ',')
            newFile.write(str(atm.atomID) + ',')
            newFile.write(str(atm.charge) + ',')
            newFile.write(str(atm.pd) + ',')
            newFile.write(str(atm.avrg_bf) + ',')
            newFile.write(str(atm.bd) + '\n')

        newFile.close()

    def make_histogram(self, threshold, highlightAtoms):  # Need to increase speed
        # Returns a kernel density plot of the B_damage values of every atom
        # considered for B_damage analysis. The area under the kernel density
        # plot (which is equal to 1) is then calculated, and a boundary is drawn
        # at the B_damage value above which a particular % (equal to the value
        # of the 'threshold' argument as defined in INPUT.txt, default value is
        # 2%) of the area under the curve lies. Those atoms which lie above this
        # threshold are listed in an html file. Any atom numbers listed in
        # highlightAtoms argument as defined in INPUT.txt will be marked on the
        # kernel density plot.

        import matplotlib.pyplot as plt
        import seaborn as sns
        import pandas as pd

        # Generates kernel density plot
        plt.clf()  # Prevents kernel density plots of all atoms considered for
        # B_damage analysis, and the subset of atoms considered for calculation
        # of the global B_damage metric, from being plotted on the same axes.
        line1 = sns.distplot(self.df.BDAM.values, hist=False, rug=True)

        # Extracts an array of (x, y) coordinate pairs evenly spaced along the
        # x(B_damage)-axis. These coordinate pairs are used to calculate the
        # area under the curve via the trapezium rule. Concomitantly, the (x, y)
        # coordinate pair values are separated depending upon whether they lie
        # below or above the threshold area defined in INPUT.txt or not. The
        # threshold B_damage boundary is then taken as the smallest x coordinate
        # of the (x, y) pairs which lie above the threshold area.
        xy_values = line1.get_lines()[0].get_data()
        x_values = xy_values[0]
        y_values = xy_values[1]

        total_area = 0
        for index, value in enumerate(y_values):
            if index != len(y_values) - 1:
                area = (((y_values[int(index)] + y_values[int((index)+1)]) / 2)
                        * (float(x_values[len(x_values)-1] - x_values[0]) / float(len(x_values)-1)))
                total_area = total_area + area

        area_LHS = 0
        atoms_list = []
        for index, value in enumerate(y_values):
            if area_LHS <= (1-float(threshold)) * total_area:
                atoms_list.append(value)
                area = (((y_values[int(index)] + y_values[int((index)+1)]) / 2)
                        * (float(x_values[len(x_values)-1] - x_values[0]) / float(len(x_values)-1)))
                if area_LHS + area <= (1-float(threshold)) * total_area:
                    area_LHS = area_LHS + area
                else:
                    break

        x_values_RHS = x_values[len(atoms_list):]
        RHS_Bdam_values = []
        for value in self.df.BDAM.values:
            if value >= x_values_RHS[0]:
                RHS_Bdam_values.append(value)

        # Those atoms considered for B_damage analysis with an associated
        # B_damage value greater than the threshold boundary are listed in an
        # html file.
        df_ordered = self.df.sort_values(by='BDAM', ascending=False)
        df_trunc = df_ordered.head(len(RHS_Bdam_values))
        decimals = pd.Series([2, 2, 2], index=['BFAC', 'AVRG BF', 'BDAM'])
        df_trunc = df_trunc.round(decimals)
        df_trunc.to_html(str(self.pdb_file_path) + '_Bdamage.html', index=False,
                         float_format='%11.3f')

        # Marks the position of the threshold boundary, plus the positions of any
        # atoms whose numbers are listed in highlightAtoms argument as defined in
        # INPUT.txt, on the kernel density plot.
        highlighted_atoms = [None]

        boundary_line = plt.plot([x_values_RHS[0], x_values_RHS[0]],
                                 [0, max(y_values)], linewidth=2, color='black',
                                 label=' boundary = {:.2f}\n (threshold = {:})'.format(x_values_RHS[0],
                                 threshold))
        highlighted_atoms.append(boundary_line)

        if len(highlightAtoms) != 0:
            lines = [None]
            for atm in highlightAtoms:
                for index, value in enumerate(self.df.ATMNUM.values):
                    if float(atm) == value:
                        lines[0] = index
                for line in lines:
                    for index, value in enumerate(self.df.BDAM.values):
                        if line == index:
                            m = plt.plot([value, value], [0, max(y_values)],
                                         linewidth=2, label=' atom ' + str(atm)
                                         + '\n B_damage = {:.2f}'.format(value))
                            highlighted_atoms.append(m)

        plt.legend(handles=highlighted_atoms[0])
        plt.xlabel('B Damage')
        plt.ylabel('Normalised Frequency')
        plt.title(str(self.pdb_code) + ' kernel density plot')
        plt.savefig(str(self.pdb_file_path)+"_Bdamage.png")

    def make_colourbyBdam_pdb(self, header_lines, footer_lines, atomList):
        # Writes pdb files in which B factor values are replaced by B_damage
        # values to allow structure when viewed with molecular graphics software
        # to be coloured by B_damage.

        import numpy as np

        # Writes PDB file in which every atom can be coloured by its B_damage
        # values.
        newPDBfile = open(str(self.pdb_file_path) + '_Bdamage.pdb', 'a')

        for line in header_lines:
            newPDBfile.write(line)

        for atm in atomList:
            a = str(atm.lineID)
            b = int(atm.atomNum)
            c = str(atm.atomType)
            d = str(atm.resiType)
            e = str(atm.chainID)
            f = int(atm.resiNum)
            g = float(atm.xyzCoords[0][0])
            h = float(atm.xyzCoords[1][0])
            j = float(atm.xyzCoords[2][0])
            k = float(atm.occupancy)
            l = float(np.log(atm.bd))  # Converts log normal B_damage distribution
            # to a linear distribution such that a fixed change in colour will
            # represent a fixed change in B_damage.
            m = str(atm.atomID)
            n = str(atm.charge)
            newLine = '%-6s%5d  %-3s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n' % (a, b, c, d, e, f, g, h, j, k, l, m, n)
            newPDBfile.write(newLine)

        for line in footer_lines:
            newPDBfile.write(line)

        newPDBfile.close()

    def calculate_Bnet(self, window_name, pdt_name):
        # Plots a kernel density estimate of Cys S, Glu O and Asp O atoms from
        # the subset of atoms considered for B_damage analysis. The Bnet metric
        # is then calculated as the ratio of the areas under the curve either side
        # of 1.

        import os
        import matplotlib.pyplot as plt
        import seaborn as sns
        import pandas as pd

        # Calculates median
        median = self.df.BDAM.median()

        # Selects atoms of Glu / Asp terminal oxygens from complete DataFrame.
        a = self.df[self.df.RESNAME.isin(['GLU'])][self.df.ATMNAME.isin(['OE1', 'OE2'])]
        b = self.df[self.df.RESNAME.isin(['ASP'])][self.df.ATMNAME.isin(['OD1', 'OD2'])]
        dataframes = [a, b]
        prot = pd.concat(dataframes)

        # Selects atoms of sugar-phosphate C-O bonds from complete DataFrame.
        na = self.df[self.df.ATMNAME.isin(["O3'", "O5'", "C3'", "C5'"])]

        if prot.empty and na.empty:
            print('\nNo sites used for Bnet calculation present in structure\n')
            pass

        if not prot.empty:
            plt.clf()  # Prevents kernel density plots of all atoms considered for
            # B_damage analysis, and the subset of atoms considered for calculation
            # of the global B_damage metric, from being plotted on the same axes.
            plot = sns.distplot(prot.BDAM.values, hist=False, rug=True)
            plt.xlabel("B Damage")
            plt.ylabel("Normalised Frequency")
            plt.title(str(self.pdb_code) + ' kernel density plot')

            # Extracts an array of (x, y) coordinate pairs evenly spaced along
            # the x(B_damage)-axis from the kernel density plot. These coordinate
            # pairs are used to calculate, via the trapezium rule, the area under
            # the curve between the smallest value of x and 1 (= area LHS), and
            # the area under the curve between 1 and the largest value of x
            # (= area RHS). The global B_damage metric is then calculated as the
            # ratio of area RHS to area LHS.
            xy_values = plot.get_lines()[0].get_data()
            x_values = xy_values[0]
            y_values = xy_values[1]

            # Calculates area RHS
            x_values_RHS = x_values[x_values >= median]
            x_min_index = len(x_values) - len(x_values_RHS)
            y_values_RHS = y_values[x_min_index:]

            x_max_RHS = int(len(x_values_RHS) - 1)
            x_distance_RHS = x_values_RHS[x_max_RHS] - x_values_RHS[0]

            total_area_RHS = 0

            for index, value in enumerate(y_values_RHS):
                if float(index) != (len(y_values_RHS) - 1):
                    area_RHS = (((y_values_RHS[int(index)] + y_values_RHS[int((index) + 1)]) / 2)
                                * (float(x_distance_RHS) / float(x_max_RHS)))
                    total_area_RHS = total_area_RHS + area_RHS

            # Calculates area LHS
            x_values_LHS = x_values[x_values <= median]
            x_max_index = len(x_values_LHS) - 1
            y_values_LHS = y_values[:x_max_index]

            x_max_LHS = int(len(x_values_LHS) - 1)
            x_distance_LHS = x_values_LHS[x_max_index] - x_values_LHS[0]

            total_area_LHS = 0

            for index, value in enumerate(y_values_LHS):
                if float(index) != (len(y_values_LHS) - 1):
                    area_LHS = (((y_values_LHS[int(index)] + y_values_LHS[int((index) + 1)]) / 2)
                                * (float(x_distance_LHS) / float(x_max_LHS)))
                    total_area_LHS = total_area_LHS + area_LHS

            # Calculates area ratio ( = global B_damage metric)
            ratio = total_area_RHS / total_area_LHS

            plt.annotate('Bnet = {:.2f}'.format(ratio),
                         xy=((max((plot.get_lines()[0].get_data())[0])*0.65),
                         (max(y_values)*0.9)))  # Need to fix annotation position
            plt.annotate('Median = {:.2f}'.format(median),
                         xy=((max((plot.get_lines()[0].get_data())[0])*0.65),
                         (max(y_values)*0.8)))  # Need to fix annotation position
            plt.savefig(str(self.pdb_file_path)+'_Bnet_Protein.png')

            if not os.path.isfile('Logfiles/Bnet_Protein.csv'):
                Bnet_list = open('Logfiles/Bnet_Protein.csv', 'w')
                Bnet_list.write('PDB' + ',')
                Bnet_list.write('Bnet' + ',')
                Bnet_list.write('Window_size (%)' + ',')
                Bnet_list.write('PDT' + ',')
                Bnet_list.close()
            Bnet_list = open('Logfiles/Bnet_Protein.csv', 'a')
            Bnet_list.write('\n%s' % self.pdb_code + ',')
            Bnet_list.write('%s' % ratio + ',')
            Bnet_list.write('%s' % window_name + ',')
            Bnet_list.write('%s' % pdt_name + ',')
            Bnet_list.close()

        if not na.empty:
            plt.clf()  # Prevents kernel density plots of all atoms considered for
            # B_damage analysis, and the subset of atoms considered for calculation
            # of the global B_damage metric, from being plotted on the same axes.
            plot = sns.distplot(na.BDAM.values, hist=False, rug=True)
            plt.xlabel("B Damage")
            plt.ylabel("Normalised Frequency")
            plt.title(str(self.pdb_code) + ' kernel density plot')

            # Extracts an array of (x, y) coordinate pairs evenly spaced along
            # the x(B_damage)-axis from the kernel density plot. These coordinate
            # pairs are used to calculate, via the trapezium rule, the area under
            # the curve between the smallest value of x and 1 (= area LHS), and
            # the area under the curve between 1 and the largest value of x
            # (= area RHS). The global B_damage metric is then calculated as the
            # ratio of area RHS to area LHS.
            xy_values = plot.get_lines()[0].get_data()
            x_values = xy_values[0]
            y_values = xy_values[1]

            # Calculates area RHS
            x_values_RHS = x_values[x_values >= median]
            x_min_index = len(x_values) - len(x_values_RHS)
            y_values_RHS = y_values[x_min_index:]

            x_max_RHS = int(len(x_values_RHS) - 1)
            x_distance_RHS = x_values_RHS[x_max_RHS] - x_values_RHS[0]

            total_area_RHS = 0

            for index, value in enumerate(y_values_RHS):
                if float(index) != (len(y_values_RHS) - 1):
                    area_RHS = (((y_values_RHS[int(index)] + y_values_RHS[int((index) + 1)]) / 2)
                                * (float(x_distance_RHS) / float(x_max_RHS)))
                    total_area_RHS = total_area_RHS + area_RHS

            # Calculates area LHS
            x_values_LHS = x_values[x_values <= median]
            x_max_index = len(x_values_LHS) - 1
            y_values_LHS = y_values[:x_max_index]

            x_max_LHS = int(len(x_values_LHS) - 1)
            x_distance_LHS = x_values_LHS[x_max_index] - x_values_LHS[0]

            total_area_LHS = 0

            for index, value in enumerate(y_values_LHS):
                if float(index) != (len(y_values_LHS) - 1):
                    area_LHS = (((y_values_LHS[int(index)] + y_values_LHS[int((index) + 1)]) / 2)
                                * (float(x_distance_LHS) / float(x_max_LHS)))
                    total_area_LHS = total_area_LHS + area_LHS

            # Calculates area ratio ( = global B_damage metric)
            ratio = total_area_RHS / total_area_LHS

            plt.annotate('Bnet = {:.2f}'.format(ratio),
                         xy=((max((plot.get_lines()[0].get_data())[0])*0.65),
                         (max(y_values)*0.9)))  # Need to fix annotation position
            plt.annotate('Median = {:.2f}'.format(median),
                         xy=((max((plot.get_lines()[0].get_data())[0])*0.65),
                         (max(y_values)*0.8)))  # Need to fix annotation position
            plt.savefig(str(self.pdb_file_path)+'_Bnet_NA.png')

            if not os.path.isfile('Logfiles/Bnet_NA.csv'):
                Bnet_list = open('Logfiles/Bnet_NA.csv', 'w')
                Bnet_list.write('PDB' + ',')
                Bnet_list.write('Bnet' + ',')
                Bnet_list.write('Window_size (%)' + ',')
                Bnet_list.write('PDT' + ',')
                Bnet_list.close()
            Bnet_list = open('Logfiles/Bnet_NA.csv', 'a')
            Bnet_list.write('\n%s' % self.pdb_code + ',')
            Bnet_list.write('%s' % ratio + ',')
            Bnet_list.write('%s' % window_name + ',')
            Bnet_list.write('%s' % pdt_name + ',')
            Bnet_list.close()
