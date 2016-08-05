

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
        #concatenate temporary variables into a single string with correct PDB formatting
        newLine = '%-6s%5d  %-3s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n' % (a,b,c,d,e,f,g,h,j,k,l,m,n)
        #write line to new PDB file
        newPDBfile.write(newLine)
    for line in eof:
        # write line to new PDB file
        newPDBfile.write(line)
    print 'New PDB file saved to %s' % newPDBfilename
    newPDBfile.close()
# end makePDB


def writeBdam(atomList, noAtm, avB, binSize, minPD, adjNo, bdamatomList):
    # method to write an output file for the calculated Bdamage values
    import pandas as pd
    import numpy as np

    numberOfRows = len(bdamatomList)
    df = pd.DataFrame(index=np.arange(0, numberOfRows),
                      columns=['REC',
                               'ATMNUM',
                               'ATMNAME',
                               'RESNAME',
                               'CHAIN',
                               'RESNUM',
                               'XPOS',
                               'YPOS',
                               'ZPOS',
                               'OCC',
                               'BFAC',
                               'ELEMENT',
                               'CHARGE',
                               'PD',
                               'BIN',
                               'GNUM',
                               'ANUM',
                               'AVRG BF',
                               'BDAM'])

    for atm in atomList:
        # take object information to a set of temporary variables

        q = int(atm.gn)
        gNo = q - 1
        binMin = int(adjNo + gNo*binSize)
        binMax = int(adjNo + q*binSize)

        dfatm = pd.DataFrame([[atm.lineID,
                               atm.atomNum,
                               atm.atomType,
                               atm.resiType,
                               atm.chainID,
                               atm.resiNum,
                               atm.xyzCoords[0][0],
                               atm.xyzCoords[1][0],
                               atm.xyzCoords[2][0],
                               atm.occupancy,
                               atm.bFactor,
                               atm.atomID,
                               atm.charge,
                               atm.pd,
                               ' %3d <= PD < %-3d' % (binMin, binMax),
                               atm.gn,
                               noAtm[gNo],
                               avB[gNo],
                               atm.bd]],
                             columns=['REC',
                                      'ATMNUM',
                                      'ATMNAME',
                                      'RESNAME',
                                      'CHAIN',
                                      'RESNUM',
                                      'XPOS',
                                      'YPOS',
                                      'ZPOS',
                                      'OCC',
                                      'BFAC',
                                      'ELEMENT',
                                      'CHARGE',
                                      'PD',
                                      'BIN',
                                      'GNUM',
                                      'ANUM',
                                      'AVRG BF',
                                      'BDAM'])

        df = df.append(dfatm, ignore_index=True)

    return df
