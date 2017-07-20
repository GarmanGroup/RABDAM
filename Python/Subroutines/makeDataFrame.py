

def makePDB(bof, atomList, eof, newPDBfilename):
    # Writes a pdb file containing a complete set of atom information for all
    # atoms in 'atomList', plus header and footer information.

    newPDBfile = open(newPDBfilename, 'a')

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
        l = float(atm.bFactor)
        m = str(atm.atomID)
        n = str(atm.charge)
        # Atom properties are appropriately ordered and spaced, and reported
        # to the expected number of significant figures, for the PDB file
        # format.
        newLine = '%-6s%5d  %-3s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n' % (
            a, b, c, d, e, f, g, h, j, k, l, m, n
            )
        newPDBfile.write(newLine)

    for line in eof:
        newPDBfile.write(line)

    print '\nNew PDB file saved to %s' % newPDBfilename
    newPDBfile.close()


def writeDataFrame(bdamatomList):
    # Returns a DataFrame containing a complete set of atom information
    # (including both that provided in the input PDB file and also the
    # B_damage values calculated by RABDAM) for all atoms considered for
    # B_damage analysis.

    import pandas as pd

    # Initialises a list for each atom property considered.
    REC = [None]*len(bdamatomList)
    ATMNUM = [None]*len(bdamatomList)
    ATMNAME = [None]*len(bdamatomList)
    RESNAME = [None]*len(bdamatomList)
    CHAIN = [None]*len(bdamatomList)
    RESNUM = [None]*len(bdamatomList)
    XPOS = [None]*len(bdamatomList)
    YPOS = [None]*len(bdamatomList)
    ZPOS = [None]*len(bdamatomList)
    OCC = [None]*len(bdamatomList)
    BFAC = [None]*len(bdamatomList)
    ELEMENT = [None]*len(bdamatomList)
    CHARGE = [None]*len(bdamatomList)
    PD = [None]*len(bdamatomList)
    AVRG_BF = [None]*len(bdamatomList)
    BDAM = [None]*len(bdamatomList)

    # Lists are filled with the relevant values of the properties associated
    # with each of the atoms considered for B_damage analysis.
    for index, atm in enumerate(bdamatomList):
        REC[index] = atm.lineID
        ATMNUM[index] = atm.atomNum
        ATMNAME[index] = atm.atomType
        RESNAME[index] = atm.resiType
        CHAIN[index] = atm.chainID
        RESNUM[index] = atm.resiNum
        XPOS[index] = atm.xyzCoords[0][0]
        YPOS[index] = atm.xyzCoords[1][0]
        ZPOS[index] = atm.xyzCoords[2][0]
        OCC[index] = atm.occupancy
        BFAC[index] = atm.bFactor
        ELEMENT[index] = atm.atomID
        CHARGE[index] = atm.charge
        PD[index] = atm.pd
        AVRG_BF[index] = atm.avrg_bf
        BDAM[index] = atm.bd

    # Lists are concatenated into the colummns of a DataFrame.
    df = pd.DataFrame({'REC': REC,
                       'ATMNUM': ATMNUM,
                       'ATMNAME': ATMNAME,
                       'RESNAME': RESNAME,
                       'CHAIN': CHAIN,
                       'RESNUM': RESNUM,
                       'XPOS': XPOS,
                       'YPOS': YPOS,
                       'ZPOS': ZPOS,
                       'OCC': OCC,
                       'BFAC': BFAC,
                       'ELEMENT': ELEMENT,
                       'CHARGE': CHARGE,
                       'PD': PD,
                       'AVRG BF': AVRG_BF,
                       'BDAM': BDAM})

    # DataFrame columns are ordered.
    cols = df.columns.tolist()
    cols_a = [cols[10]] + [cols[1]] + [cols[0]] + [cols[11]] + [cols[5]]
    cols_b = [cols[12]] + [cols[13]] + [cols[14]] + [cols[15]] + [cols[8]]
    cols_c = [cols[4]] + [cols[7]] + [cols[6]] + [cols[9]]
    cols_d = [cols[2]] + [cols[3]]
    cols = cols_a + cols_b + cols_c + cols_d
    df = df[cols]

    return df
