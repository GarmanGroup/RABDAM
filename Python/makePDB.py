# Copyright Thomas Dixon 2015
def makePDB(atomList, newPDBfilename):
#    import os
    from parsePDB import atom
#    if os.path.exists(newPDBfilename):
#        print 'File %s already exists' % newPDBfilename
#        return
    newPDBfile = open(newPDBfilename,'a')
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
        newLine = '%-6s%5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s' % (a,b,c,d,e,f,g,h,j,k,l,m,n)
        newPDBfile.write(newLine)
    newPDBfile.close()
#end makePDB