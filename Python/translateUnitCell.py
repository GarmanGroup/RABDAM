#Copyright Thomas Dixon 2015
#convert the unit cell paramenters to Cartesian coordinates
def convertToCartesian(unitCell):
    import math #facilitates using more intricate maths operations
    import numpy as np #facilitates matrix manipulation
    print 'Converting the unit cell basis vectors to Cartesian coordinates'
    a = float(unitCell[0])
    b = float(unitCell[1])
    c = float(unitCell[2])
    alpha = float(unitCell[3])
    beta = float(unitCell[4])
    gamma = float(unitCell[5])
    #define parameter v, which is the volume of a unit parallelepiped with the 
    #same angles as the unit cell
    v = math.sqrt(1 - math.pow((math.cos(alpha)),2) - math.pow((math.cos(beta)),2) - math.pow((math.cos(gamma)),2) + 2*(math.cos(alpha))*(math.cos(beta))*(math.cos(gamma)))
    #define the elements of the conversion matrix
    a11 = a
    a12 = b*math.cos(gamma)
    a13 = c*math.cos(beta)
    a21 = 0
    a22 = b*math.sin(gamma)
    a23 = c*(math.cos(alpha)-(math.cos(beta)*math.cos(gamma))/math.sin(gamma))
    a31 = 0
    a32 = 0
    a33 = c*(v/math.sin(gamma))
    #create the conversion matrix
    conversionMatrix = np.array([[a11,a12,a13],[a21,a22,a23],[a31,a32,a33]])
    #define the fractional basis vectors for each direction; a, b and c
    aVector = np.array([[1],[0],[0]])
    bVector = np.array([[0],[1],[0]])
    cVector = np.array([[0],[0],[1]])
    #convert the lattice basis vector in each direction to Cartesian coordinates
    aCartesianVector = np.dot(conversionMatrix, aVector)
    bCartesianVector = np.dot(conversionMatrix, bVector)
    cCartesianVector = np.dot(conversionMatrix, cVector)
    cartesianVectors = ([aCartesianVector],[bCartesianVector],[cCartesianVector])
    return cartesianVectors
    print 'Conversion complete\n'
#end translateUnitCell
    
def translateUnitCell(atomList, cartesianVectors, atrans, btrans, ctrans):
    from parsePDB import atom #for utilising the 'atom' class
    import numpy as np #facilitates matrix manipulation
    #create puppet list to fill with atom objects
    transAtoms = [] 
    #loop through all atoms in supplied list
    for atm in xrange (0,len(atomList)):
        #obtain atom information for atom n from atomList
        n = atom(atomList(atm))
        #turn the xyzCoords of atom n into a matrix
        cartCoords = np.array(n.xyzCoords)
        #create a matrix of the number of unit cells to translate in a,b,c axes
        fracVectorTrans = np.array([atrans, btrans, ctrans])
        #multiply the above matrices to give a Cartesian translation
        vectorTrans = np.dot(fracVectorTrans, cartesianVectors)
        #apply this transformation to the atoms xyzCoords and write back to atom n
        n.xyzCoords = np.dot(vectorTrans, cartCoords)
        #append the translated atom object to list
        transAtoms.append(n)
    return transAtoms   
#end TranslateUnitCell