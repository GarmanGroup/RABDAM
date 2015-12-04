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
        g = (atm.xyzCoords)
        h = g[1]
        h = h[0]
        j = g[2]
        j = j[0]
        g = g[0]
        g = g[0]
        k = float(atm.occupancy)
        l = float(atm.bFactor)
        m = str(atm.atomID)
        n = str(atm.charge)
        #concatenate temporary variables into a single string with correct PDB formatting
        newLine = '%-6s%5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n' % (a,b,c,d,e,f,g,h,j,k,l,m,n)
        #write line to new PDB file
        newPDBfile.write(newLine)
    for line in eof:
        #write line to new PDB file
        newPDBfile.write(line)
    print 'New PDB file saved to %s' % newPDBfilename
    newPDBfile.close()
#end makePDB