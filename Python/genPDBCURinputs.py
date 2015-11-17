#Copyright Thomas Dixon 2015
def genPDBCURinputs(pathToPDB, PDBdirectory, delhydrogen):
    #import os for OS usability
    import os
    #create path to PDBCUR input file
    PDBCURinputFile = '%s/PDBCURinput.txt' % PDBdirectory 
    #check if an input file has already been created
    if os.path.exist(PDBCURinputFile):
        #inform user file already exists
        print 'Input file for PDBCUR already exist at %s/PDBCURinputFile.txt' % PDBdirectory
        return
        
    #open a text file for writing
    
    
    
#end