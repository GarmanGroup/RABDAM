function atomicCoordinates = translateUnitCellAlternative(atomicCoordinates,xyzAtom,aBasis,bBasis,cBasis,translationVector,numberOfAtoms)

%Calculate the total translation vector
totalTranslationVector = aBasis*translationVector(1) + bBasis*translationVector(2) + cBasis*translationVector(3);

%Create array array same size as xyzAtom that contains the total
%translation vector in every column
totalTranslationVector2dArray = repmat(totalTranslationVector',numberOfAtoms,1);

%Translate each atom.
xyzTranslatedAtom = xyzAtom + totalTranslationVector2dArray;

%Replace the atom xyz coordinates to the translated atom coordinates
atomicCoordinates(:,9:11) = num2cell(xyzTranslatedAtom);

end
