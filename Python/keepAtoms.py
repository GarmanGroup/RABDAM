def keepAtoms(atomList, keepHETATM, keepHOH)
	from parsePDB import atom as a
	if keepHETATM = True
		if keepHOH = True
			return atomList
	for atom in atomList:
		if keepHETATM = False:
			if atom.atomType = 'HETATM'
				atom.pop
		elif keepHOH = False:
			if atom.resiType = 'HOH'
				atom.pop
			elif atom.resiType = 'WAT'
				atom.pop
			elif atom.resiType = 'H2O'
				atom.pop
			elif atom.resiType = 'DOD'
				atom.pop
	return atomList
#end keepAtoms