

# method to make a pdb from complete set of atom information and pdb header/footer
def makePDB(bof, atomList, eof, newPDBfilename):
    newPDBfile = open(newPDBfilename, 'a')
    for line in bof:
        # write line to new PDB file
        newPDBfile.write(line)
    for atm in atomList:
        # take object information to a set of temporary variables
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
        # concatenate temporary variables into a single string with correct PDB formatting
        newLine = '%-6s%5d  %-3s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n' % (a, b, c, d, e, f, g, h, j, k, l, m, n)
        # write line to new PDB file
        newPDBfile.write(newLine)
    for line in eof:
        # write line to new PDB file
        newPDBfile.write(line)
    print 'New PDB file saved to %s' % newPDBfilename
    newPDBfile.close()
# end makePDB


def writeBdam(noAtm, avB, binSize, minPD, adjNo, bdamatomList):
    # method to write calculated Bdamage values to a dataframe
    import pandas as pd

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
    BIN = [None]*len(bdamatomList)
    GNUM = [None]*len(bdamatomList)
    ANUM = [None]*len(bdamatomList)
    AVRG_BF = [None]*len(bdamatomList)
    BDAM = [None]*len(bdamatomList)

    for index, atm in enumerate(bdamatomList):
        q = int(atm.gn)
        gNo = q - 1
        binMin = int(adjNo + gNo*binSize)
        binMax = int(adjNo + q*binSize)

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
        BIN[index] = '%3d <= PD < %-3d' % (binMin, binMax)
        GNUM[index] = atm.gn
        ANUM[index] = noAtm[gNo]
        AVRG_BF[index] = avB[gNo]
        BDAM[index] = atm.bd

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
                       'BIN': BIN,
                       'GNUM': GNUM,
                       'ANUM': ANUM,
                       'AVRG BF': AVRG_BF,
                       'BDAM': BDAM})
    cols = df.columns.tolist()
    cols_a = [cols[13]] + [cols[2]] + [cols[1]] + [cols[14]] + [cols[7]]
    cols_b = [cols[15]] + [cols[16]] + [cols[17]] + [cols[18]] + [cols[11]]
    cols_c = [cols[5]] + [cols[9]] + [cols[8]] + [cols[12]] + [cols[6]]
    cols_d = [cols[10]] + [cols[0]] + [cols[3]] + [cols[4]]
    cols = cols_a + cols_b + cols_c + cols_d
    df = df[cols]
    return df
