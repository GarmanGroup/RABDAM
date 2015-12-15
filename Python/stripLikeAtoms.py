def stripLikeAtoms(atomList1,atomList2)
	from parsePDB import atom as a
	likeAtomList = []
	for atom in atomList1,
		for atm in atomList2
			if atom.ID = atm.ID
			   && atom.resType = atom.resType
			   && ...
				likeAtomList.append(atom)
				atom.pop
				atm.pop
			    break
		continue
	return likeAtomList
#end stripLikeAtoms