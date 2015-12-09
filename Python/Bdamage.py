# Copyright Thomas Dixon 2015

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

#Calculate packing density for all atoms in the original PDB file
def calcPackingDensity(AUatomList, atomList, PDT):
    from parsePDB import atom as a #for utilising the 'atom' class
    from Bdamage import calcDist #calcualtes the distance between two atoms in 3D space
    #for every atom in the asymmetric unit
    for atmAU in AUatomList:
        #set packing density counter to 0
        PD = 0
        #for every retained atom
        for atm in atomList:
            #calcualte the distance between the two atoms
            if calcDist(atmAU, atm) < PDT:
                #if the distance is less than the PDT, increment the counter
                PD = PD + 1
        #return 
        atmAU.PD = PD
    return AUatomList
#end calcPackingDensity