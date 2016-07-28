# write graphs


def make_csv(atomList, filename, noAtm, avB, binSize, minPD, adjNo):
    import os
    if os.path.exists(filename):
        return
    else:
        newFile = open(filename, 'w')

        newFile.write('REC      = RECORD NAME\n'
                      'ATMNUM   = ATOM SERIAL NUMBER\n'
                      'ATMNAME  = ATOM NAME\n'
                      'RESNAME  = RESIDUE NAME\n'
                      'CHAIN    = CHAIN IDENTIFIER\n'
                      'RESNUM   = RESIDUE SEQUENCE NUMBER\n'
                      'XPOS     = ORTHOGONAL COORDINATES FOR X IN ANGSTROMS\n'
                      'YPOS     = ORTHOGONAL COORDINATES FOR Y IN ANGSTROMS\n'
                      'ZPOS     = ORTHOGONAL COORDINATES FOR Z IN ANGSTROMS\n'
                      'OCC      = OCCUPANCY\n'
                      'BFAC     = B FACTOR (TEMPERATURE FACTOR)\n'
                      'ELEMENT  = ELEMENT SYMBOL\n'
                      'CHARGE   = CHARGE ON ATOM\n'
                      'PD       = PACKING DENSITY (ATOMIC CONTACT NUMBER)\n'
                      'BIN      = SIMILAR PACKING DENSITY BIN\n'
                      'GNUM     = SIMILAR PACKING DENSITY ENVIRONMENT GROUP NUMBER\n'
                      'ANUM     = NUMBER OF ATOMS IN SIMILAR PACKING DENSITY GROUP\n'
                      'AVRG BF  = AVERAGE B FACTOR FOR ATOMS IN SIMILAR PACKING DENSITY ENVIRONMENT\n'
                      'BDAM     = BDAMAGE VALUE\n'
                      '\n')

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
                      'BIN' + ','
                      'GNUM' + ','
                      'ANUM' + ','
                      'AVRG BF' + ','
                      'BDAM' + '\n')
        for atm in atomList:
            q = int(atm.gn)
            gNo = q - 1
            binMin = int(adjNo + gNo*binSize)
            binMax = int(adjNo + q*binSize)

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
            newFile.write(str(' %3d <= PD < %-3d' % (binMin, binMax)) + ',')
            newFile.write(str(atm.gn) + ',')
            newFile.write(str(noAtm[gNo]) + ',')
            newFile.write(str(avB[gNo]) + ',')
            newFile.write(str(atm.bd) + '\n')

        print 'New PDB file saved to %s' % filename
        newFile.close()


def make_histogram(df, fileName, PDBcode, threshold):
    import matplotlib.pyplot as plt
    import numpy as np
    import seaborn as sns

    # generate B_damage kernel density plot
    line1 = sns.distplot(df.BDAM.values, hist=False, rug=True)

    # extract data from kernel density plot, calculate area under curve
    xy_values = line1.get_lines()[0].get_data()

    x_values = xy_values[0]
    y_values = xy_values[1]

    total_area = 0
    for index, value in enumerate(y_values):
        if index != len(y_values) - 1:
            area = (((y_values[int(index)] + y_values[int((index) + 1)]) / 2) * (float(x_values[len(x_values) - 1] - x_values[0]) / float(len(x_values) - 1)))
            total_area = total_area + area

    area_LHS = 0
    atoms_list = []
    for index, value in enumerate(y_values):
        if area_LHS <= (1 - float(threshold)) * total_area:
            atoms_list.append(value)
            area = (((y_values[int(index)] + y_values[int((index) + 1)]) / 2) * (float(x_values[len(x_values) - 1] - x_values[0]) / float(len(x_values) - 1)))
            if area_LHS + area <= (1 - float(threshold)) * total_area:
                area_LHS = area_LHS + area
            else:
                break

    x_values_RHS = x_values[len(atoms_list):]

    # write dataframe to html, highlighting atoms with B_damage values which lie above the 5% threshold
    RHS_Bdam_values = []
    for index, value in enumerate(df.BDAM.values):
        if value > x_values_RHS[0]:
            RHS_Bdam_values.append(value)

    x = df.sort_values(by='BDAM', ascending=0)
    y = x.head(len(RHS_Bdam_values))
    y.to_html(str(fileName) + 'Bdamage.html', index=False)

    # draw a line on kernel density plot at 5% threshold
    plt.plot([x_values_RHS[0], x_values_RHS[0]], [0, max(y_values)], linewidth=2)
    plt.annotate(' boundary = {:.2f}'.format(x_values_RHS[0]), xy=[x_values_RHS[0], (0.8*max(y_values))])
    plt.xlabel('B Damage')
    plt.ylabel('Frequency')
    plt.title(str(PDBcode) + ' kernel density plot')
    plt.savefig(str(fileName)+"Bdamage.png")

    return x_values_RHS


def make_colourbyBdam_pdb(df, bof, eof, fileName, atomList, x_values_RHS):
    import numpy as np

    # writes PDB file with B damage on a linear scale
    newPDBfile = open(str(fileName) + 'Bdamage.pdb', 'a')

    for line in bof:
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
        l = float(np.log(atm.bd))
        m = str(atm.atomID)
        n = str(atm.charge)
        newLine = '%-6s%5d  %-3s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n' % (a,b,c,d,e,f,g,h,j,k,l,m,n)
        newPDBfile.write(newLine)

    for line in eof:
        newPDBfile.write(line)

    newPDBfile.close()

    # writes PDB file highlighting only those atoms with B damage values above 5% threshold
    newPDBfile = open(str(fileName) + 'Bdamage_above_boundary.pdb', 'a')

    for line in bof:
        newPDBfile.write(line)

    for atm in atomList:
        if float(atm.bd) > x_values_RHS[0]:
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
            l = float(atm.bd)
            m = str(atm.atomID)
            n = str(atm.charge)
            newLine = '%-6s%5d  %-3s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n' % (a,b,c,d,e,f,g,h,j,k,l,m,n)
            newPDBfile.write(newLine)
        else:
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
            l = 0
            m = str(atm.atomID)
            n = str(atm.charge)
            newLine = '%-6s%5d  %-3s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n' % (a,b,c,d,e,f,g,h,j,k,l,m,n)
            newPDBfile.write(newLine)

    for line in eof:
        newPDBfile.write(line)

    newPDBfile.close()


def calculate_global_BDam(df):
    pass
