# Copyright 2015 Thomas Dixon
# With thanks to Jonathan Brooks-Bartlett, Charles Bury, Markus Gerstel and Elspeth Garman

#Script to calculate B-damage for protein atoms
import sys
sys.path.insert(0,'./Subroutines')
from CalculateBdamage import cambda

#read the filename from the CL input
fileCont = open(sys.argv[1],'r')
#write the string of the function arguments
functionArgs = fileCont.read()
#remove all whitespace from the inputs
functionArgs = functionArgs.replace(" ","")
functionArgs = functionArgs.replace("\n","")
functionArgs = functionArgs.replace("\r","")
#split the input file into its components
splitArgs = functionArgs.split(",")
#parse the argument string and assign values to CaMBDa inputs
pathToPDB = splitArgs[0]
#initialise the argument defaults
pdtVal = int(14)
binVal = int(10)
aucVal = False
taVal = False

for x in xrange(1, len(splitArgs)):
    if splitArgs[x][0:3] == "PDT":
        pdtArg = splitArgs[x].split("=")
        pdtVal = float(pdtArg[len(pdtArg)-1])
        if pdtVal == int(pdtVal):
            pdtVal = int(pdtVal)
    elif splitArgs[x][0:7] == "binSize":
        binArg = splitArgs[x].split("=")
        binVal = int(binArg[len(binArg)-1])
    elif splitArgs[x][0:12] == 'createAUCpdb':
        aucArg = splitArgs[x].split("=")
        aucVal = aucArg[len(aucArg)-1]
        if aucVal == 'True':
            aucVal = True
        elif aucVal == 'False':
            aucVal = False
        else:
            sys.exit("Error 00: Input file is formatted incorrectly\nRead the handbook and amend the input file\n(Try looking at the createAUCpdb argument)")
    elif splitArgs[x][0:11] == 'createTApdb':
        taArg = splitArgs[x].split("=")
        taVal = taArg[(len(taArg))-1]
        if taVal == 'True':
            taVal = True
        elif taVal == 'False':
            taVal = False
        else:
            sys.exit("Error 00: Input file is formatted incorrectly\nRead the handbook and amend the input file\n(Try looking at the createTApdb argument)")
    else:
        sys.exit("Error 00: Input file is formatted incorrectly\nRead the handbook and amend the input file")

#run fucntion with arguments from input file
cambda(pathToPDB, PDT=pdtVal, binSize=binVal, createAUCpdb=aucVal, createTApdb=taVal)