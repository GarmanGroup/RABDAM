#remove atoms from a list that lie outside of a set of given spatial parameters
def trimAtoms(atomList, params):
    from parsePDB import atom as a #for utilising the 'atom' class
    from atomCheck import isInXYZparams
    print 'Excluding the atoms that lie outside of the box'
    totalAtm = len(atomList)
    atmIndex = 0
    while atmIndex < totalAtm:
        #extract xyz coordinates from atom information
        atomXYZ = atomList[atmIndex].xyzCoords
        #if the newly considered atoms coordinates lie within the params, retain this atom
        if isInXYZparams(atomXYZ, params):
            #advance the index by 1 to red the next line
            atmIndex = atmIndex + 1
        #otherwise, discard the atom 
        else:
            #remove atom information from atomList
            atomList.pop(atmIndex)
            #reduce the number of total atoms in atomList by 1
            totalAtm = totalAtm - 1
    print '%.0f atoms have been retained\n' % totalAtm
    #output a list of retained atom objects
    return atomList
#end trimAtoms
