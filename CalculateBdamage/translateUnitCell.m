function translatedAtomicCoordinates = translateUnitCell(atomicCoordinates,aBasis,bBasis,cBasis,translationVector)

%Set the translated atomic coordinates equal to the unaltered atomic
%coordinates.
translatedAtomicCoordinates = atomicCoordinates;

%Loop over each atom in the unit cell
for eachAtom = 1 : length(atomicCoordinates)
    
    %extract the atomic Coordinates
    xyzAtom = str2double(atomicCoordinates(eachAtom,9:11))';
    
    %translate atom
    xyzTranslatedAtom = xyzAtom + aBasis*translationVector(1) + bBasis*translationVector(2) + cBasis*translationVector(3);
    
    %convert column vector to row vector and replace the coordinates in the
    %Translated Atomic Coordinates
    translatedAtomicCoordinates(eachAtom,9:11) = num2cell(xyzTranslatedAtom');
end

end

