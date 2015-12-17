# Copyright Thomas Dixon 2015

#method to make a pdb from complete set of atom information and pdb header/footer 
def makePDB(bof, atomList, eof, newPDBfilename):
    import os #for operating system usability
    from parsePDB import atom #for utilising the 'atom' class
    if os.path.exists(newPDBfilename):
        print 'File %s already exists' % newPDBfilename
        return
    newPDBfile = open(newPDBfilename,'a')
    for line in bof:
        #write line to new PDB file
        newPDBfile.write(line)
    for atm in atomList:
        #take object information to a set of temporary variables
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
        #write line to new PDB file
        newPDBfile.write(line)
    print 'New PDB file saved to %s' % newPDBfilename
    newPDBfile.close()
#end makePDB
    
#method to write an output file for the calculated Bdamage values
def writeBdam(atomList, filename, noAtm, avB):
    import os #for operating system usability
    from parsePDB import atom #for utilising the 'atom' class
    if os.path.exists(filename):
        print 'File %s already exists' % filename
        return
    newFile = open(filename,'a')
    #write preamble legend to output file
    newFile.write('REC  = RECORD NAME\n'
                  'SER  = ATOM SERIAL NUMBER\n'
                  'ATM  = ATOM NAME\n'
                  'A    = ALETERNATE LOCATION IDENTIFIER\n'
                  'RES  = RESIDUE NAME\n'
                  'C    = CHAIN IDENTIFIER\n'
                  'RS   = RESIDUE SEQUENCE NUMBER\n'
                  'IN   = CODE FOR INSERTION OF RESIDUES\n'
                  'XPOS = ORTHOGONAL COORDINATES FOR X IN ANGSTROMS\n'
                  'YPOS = ORTHOGONAL COORDINATES FOR Y IN ANGSTROMS\n'
                  'ZPOS = ORTHOGONAL COORDINATES FOR Z IN ANGSTROMS\n'
                  'OCC  = OCCUPANCY\n'
                  'BFAC = B FACTOR (TEMPERATURE FACTOR)\n'
                  'EL   = ELEMENT SYMBOL\n'
                  'CH   = CHARGE ON ATOM\n'
                  'PD   = PACKING DENSITY (ATOMIC CONTACT NUMBER)\n'
                  'BIN  = SIMILAR PACKING DENSITY BIN\n'
                  'GN   = SIMILAR PACKING DENSITY ENVIRONMENT GROUP NUMBER\n'
                  'ANUM = NUMBER OF ATOMS IN SIMILAR PACKING DENSITY GROUP\n'
                  'AVB  = AVERAGE B FACTOR FOR ATOMS IN SIMILAR PACKING DENSITY ENVIRONMENT\n'
                  'BDAM = BDAMAGE VALUE\n'
                  '\n')
    #write column headers
    newFile.write('REC     SER ATMA RES C   RES IN XPOS    YPOS    ZPOS    OCC  BFAC            EL CH PD  BIN              GN     ANUM  AV    BDAM\n')
    for atm in atomList:
        #take object information to a set of temporary variables
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
        o = int(atm.pd)
        p = str('               ')
        q = int(atm.gn)
        gNo = q - 1
        r = int(noAtm[gNo])
        s = float(avB[gNo])
        t = float(atm.bd)
        #concatenate temporary variables into a single string with correct PDB formatting
        newLine = '%-6s%5d  %-3s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s %2s %4d %15s %3d     %4d %6.2f %5.4f\n' % (a,b,c,d,e,f,g,h,j,k,l,m,n,o,p,q,r,s,t)
        #write line to new PDB file
        newFile.write(newLine)
    print 'New PDB file saved to %s' % filename
    newFile.close()
#end writeBdam