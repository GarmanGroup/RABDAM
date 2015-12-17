# Copyright Thomas Dixon 2015
# With thanks to Charles Bury

#define class of atom objects
class atom(object):
    #Initialise class for a PDB file
    def __init__(self,lineidentifier="",atomnum=0,residuenum=0,atomtype="",resitype="",
                 chainID="",xyz_coords=[],atomidentifier="",bfactor=0,occupancy=1,charge="",
                 packingdensity=0,groupnumber=0,bdamage=0):   
        self.lineID     = lineidentifier
        self.atomNum    = atomnum
        self.resiNum    = residuenum
        self.atomType   = atomtype
        self.resiType   = resitype 
        self.chainID    = chainID
        self.xyzCoords  = xyz_coords
        self.atomID     = atomidentifier
        self.bFactor    = bfactor
        self.occupancy  = occupancy
        self.charge     = charge
        self.pd         = packingdensity
        self.group      = groupnumber
        self.bd         = bdamage
    #end Initialise class
        
    #print a summary of atom info to command line 
    def getAtomSummary(self):
        summaryString = 'Chain: {}\nResidue: {}{}\nAtom type: {}'.format(self.chainID,self.resiType,self.resiNum,self.atomType)
        print summaryString
    #end getAtomSummary
#end atom class
        
#parse pdb file with name 'fileName' and return list of atoms from 'atom' class above
def parsePDB(fileName):
    import os #for operating system usability
    import sys #for terminating script when encountering errors
    from parsePDB import atom #for utilising the 'atom' class
    #check that file exists
    if not os.path.exists(fileName):
        sys.exit('Error!!\nFile name %s not found' % fileName)
    #create puppet lists to fill with atom objects
    bof = []
    atomList = [] 
    eof = []
    beforeAtoms = True
    #open correct file name for reading
    fileOpen = open(fileName,'r') 
    #read 'ATOM' and 'HETATM' lines
    for line in fileOpen.readlines():
        if (not ('ATOM  ' in str(line[0:6])) 
            or not ('HETATM' in str(line[0:6])) 
            or not ('ANISOU' in str(line[0:6])) 
            or not ('TER   ' in str(line[0:6]))):
            if beforeAtoms:
                bof.append(line)
            else:
                eof.append(line)  
        #only select information from ATOM or HETATM lines
        if ('ATOM  ' in str(line[0:6])) or ('HETATM' in str(line[0:6])):
            beforeAtoms = False
            y = atom() #make new 'atom' object here
            #get atom properties here
            y.lineID    = str(line[0:6].strip())
            y.atomNum   = int(line[6:11].strip())
            y.atomType  = str(line[12:16].strip())
            y.resiType  = str(line[17:20].strip())                       
            y.chainID   = str(line[21:22].strip())                     
            y.resiNum   = int(line[22:26].strip())  
            y.xyzCoords = [[float(line[30:38].strip())],
                           [float(line[38:46].strip())],
                           [float(line[46:54].strip())]]    
            y.occupancy = float(line[54:60].strip())                                                    
            y.bFactor   = float(line[60:66].strip())
            y.atomID    = str(line[77:79].strip())
            y.charge    = str(line[79:81].strip())
            #append new 'atom' object to list
            atomList.append(y)
    fileOpen.close() #close .pdb file after reading
    #provide feedback to user
    print 'Finished reading in atoms --> %d atoms found in' % int(len(atomList))
    print '%s' % fileName
    #return list of atoms as output of function
    return bof, atomList, eof
#end parsePDB
    
#obtain unit cell parameters from PDB file
def getUnitCellParams(fileName):
    import os #for operating system usability   
    import math #for converting degrees to radians
    #check that file exists
    if not os.path.exists(fileName):
        print 'Error!!\nFile name {} not found'.format(fileName)
        return
    #open correct file name for reading
    fileOpen = open(fileName,'r') 
    #find and read the CRYST1 line
    for line in fileOpen.readlines():
        #only select infoprmation from CRYST1 line
        if('CRYST1' in str(line[0:6])):
            params = line.split()
            #create object 'y' with attributes for each of the Unit Cell Parameters
            a = float(params[1])
            b = float(params[2])
            c = float(params[3])
            alpha = math.radians(float(params[4]))
            beta = math.radians(float(params[5]))
            gamma = math.radians(float(params[6]))
            break
    #provide feedback to user
    print 'Unit cell parameters successfully extracted'
    fileOpen.close() #close .pdb file
    return (a, b, c, alpha, beta, gamma)
#end getUnitCellParams
    
#get minimum and maximum xyz coordinates of the asymmetric unit
def getAUparams(atomList):
    from parsePDB import atom as a #for utilising the 'atom' class
    #initialise xyz minima and maxima
    xMin = float(atomList[0].xyzCoords[0][0])
    xMax = float(atomList[0].xyzCoords[0][0])
    yMin = float(atomList[0].xyzCoords[1][0])
    yMax = float(atomList[0].xyzCoords[1][0])
    zMin = float(atomList[0].xyzCoords[2][0])
    zMax = float(atomList[0].xyzCoords[2][0])
    #loop through all atoms in input list
    for atm in atomList:
        #extract xyz coordinates from atom information
        x = float(atm.xyzCoords[0][0])
        y = float(atm.xyzCoords[1][0])
        z = float(atm.xyzCoords[2][0])
        #if the newly considered atoms coordinates lie outside of the previous 
        #max/minima, replace the relevant parameter(s)
        if x < xMin:
            xMin = x
        elif x > xMax:
            xMax = x
        if y < yMin:
            yMin = y
        elif y > yMax:
            yMax = y
        if z < zMin:
            zMin = z
        elif z < zMax:
            zMax = z
    #combine all parameters into a single list
    auParams = [xMin, xMax, yMin, yMax, zMin, zMax]
    return auParams
#end getAUparams
    
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
    
def trimAtoms2(atomList, params):
    from parsePDB import atom as a #for utilising the 'atom' class
    from atomCheck import isInXYZparams
    print 'Excluding the atoms that lie outside of the box'
    trimAtomList = []
    for atom in atomList:
        atomXYZ = atom.xyzCoords
        if isInXYZparams(atomXYZ, params):
            trimAtomList.append(atom)
    print '%.0f atoms have been retained\n' % int(len(trimAtomList))
    #output a list of retained atom objects
    return trimAtomList
#end trimAtoms2
            