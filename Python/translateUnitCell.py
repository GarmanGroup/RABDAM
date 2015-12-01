#Copyright Thomas Dixon 2015
#convert the unit cell paramenters to Cartesian coordinates
def convertToCartesian(unitCell):
    import math #facilitates using more intricate maths operations
    import numpy as np #facilitates matrix manipulation
    print 'Converting the unit cell basis vectors to Cartesian coordinates'
    a = unitCell(0)
    b = unitCell(1)
    c = unitCell(2)
    alpha = unitCell(3)
    beta = unitCell(4)
    gamma = unitCell(5)
    #define parameter v, whcih is the volume fo the unit cell
    v = math.sqrt(1-(math.cos(alpha))^2-(math.cos(beta))^2-(math.cos(gamma))^2
                  +2*(math.cos(alpha))*(math.cos(beta))*(math.cos(gamma)))
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
    conversionMatrix = np.array([a11,a12,a13],[a21,a22,a23],[a31,a32,a33])
    #define the fractional basis vectors for each direction; a, b and c
    aVector = np.array([1],[0],[0])
    bVector = np.array([0],[1],[0])
    cVector = np.array([0],[0],[1])
    #convert the lattice basis vector in each direction to Cartesian coordinates
    aCartesianVector = np.multiply(conversionMatrix, aVector)
    bCartesianVector = np.multiply(conversionMatrix, bVector)
    cCartesianVector = np.multiply(conversionMatrix, cVector)
    cartesianVectors = ([aCartesianVector],[bCartesianVector],[cCartesianVector])
    return cartesianVectors
    print unitCell
#end translateUnitCell