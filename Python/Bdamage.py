# Copyright Thomas Dixon 2015

#calculate the distance between two atoms
def calcDist(atomA, atomB):
    from parsePDB import atom as a #for utilising the 'atom' class
    import math #for obtaining the square root
    #obtain two sets of atom coordinates for the atoms
    Acoords = atomA.xyzCoords
    Bcoords = atomB.xyzCoords
    #find the linear displaements of the atoms from each other along the three Cartesian axes
    x = Acoords[0][0] - Bcoords[0][0]
    y = Acoords[1][0] - Bcoords[1][0]
    z = Acoords[2][0] - Bcoords[2][0]
    #find the distance between the two atoms 
    d = math.sqrt(x**2 + y**2 + z**2)
    return d
#end calcDist

def countInRadius(atm, atomList, r):
    from parsePDB import atom as a #for utilising the 'atom' class
    from Bdamage import calcDist #calcualtes the distance between two atoms in 3D space
    #set packing density counter to 0
    PD = 0
    #for every atom
    for atm in atomList:
        #calcualte the distance between the two atoms
        if calcDist(atm, atm) < r:
            #if the distance is less than the PDT, increment the counter
            PD = PD + 1
    #return packing density of the atom once all comparisons have been made
    atm.pd = PD
    return atm.pd
        
#Calculate packing density for all atoms in the original PDB file
def calcPDT(auAtomList, atomList, PDT):
    from parsePDB import atom as a #for utilising the 'atom' class
    from Bdamage import countInRadius #counts atoms within PDT
    #set initial values for min/maxPD
    firstPD = countInRadius(auAtomList[0], atomList, PDT)
    minPD = firstPD
    maxPD = firstPD
    #for every atom in the asymmetric unit
    for atm in xrange (1,len(auAtomList)):
        auAtm = auAtomList[atm]
        auAtm.PD = countInRadius(auAtm, atomList, PDT)
        #update min/maxPD if necessary
        if auAtm.pd < minPD:
            minPD = auAtm.PD
        elif auAtm.pd > maxPD:
            maxPD = auAtm.pd
    print 'Packing Density (PD) values successfully calculated'
    return auAtomList, minPD, maxPD
#end calcPackingDensity
    
#Segregate atoms into bins based on PD
def binAtoms(atomList, binSize, minPD, maxPD):
    from parsePDB import atom as a #for utilising the 'atom' class
    import math #to utilise more intricate maths functions
    #create value for 'adjustment number' which is a factor to be taken off all
    #PDs in order to define their group number by a ceiling function
    adjtNo = (math.floor(minPD/binSize))*binSize
    #for every atom in the atom list
    for atm in atomList:
        #obtain the PD value
        actlPD = atm.pd
        #reduce PD by the adjustment value
        adjdPD = actlPD - adjtNo + 1
        #define group number as the ceiling of adjdPD divided by bin size
        groupNo = math.ceil(adjdPD/binSize)
        atm.gn = groupNo
    noOfGroups = int(math.ceil((maxPD-adjtNo)/binSize))
    return atomList, noOfGroups
#end binAtoms
    
#calculate Bdamage value for every atom in AU
def calcBdam(atomList, numberOfGroups):
    from parsePDB import atom as a #for utilising the 'atom' class
    #initialise variables for all group numbers
    sumB = [0] * numberOfGroups
    noAtm = [0 for col in range(numberOfGroups)]
    avB = [0 for col in range(numberOfGroups)]
    #find sum of all B factors of atoms in their groups
    for atom in atomList:
        gNo = int(atom.gn - 1) #take away 1 to account for cardinality vs ordinality
        sumB[gNo] = float(sumB[gNo]) + float(atom.bFactor)
        noAtm[gNo] = int(noAtm[gNo]) + 1
    #find the average B factor for each group number
    for gNo in xrange(numberOfGroups-1):
        avB[gNo] = float(sumB[gNo])/int(noAtm[gNo])
    #calculate B damage for each atom and update this value for the atom object
    for atom in atomList:
        gNo = int(atom.gn - 1)
        atom.bd = float(atom.bFactor)/float(avB[gNo])
    #return outputs of the script
    return atomList, noAtm, avB
#end calcBdam        