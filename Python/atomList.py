# Copyright Thomas Dixon 2015

def remAtoms(atomList, remHETATM, remWaters):
    from parsePDB import atom as a #for utilising the 'atom' class
    #if the function is not supposed to remove any atoms, return unaltered atomList
    if not remHETATM:
        if not remWaters:
            return atomList
    #loop through all atoms in input list
    totalAtm = len(atomList)
    atmIndex= 0
    keptAtoms = 0
    while atmIndex < totalAtm:
        if remHETATM:
            if atom.atomType == 'HETATM':
                atomList.pop(atmIndex)
                totalAtm = totalAtm - 1
                break
            else:
                atmIndex = atmIndex + 1
                keptAtoms = keptAtoms + 1
        if remWaters:
            
#end remAtoms    