%% Function to remove atoms outside of box
function newAtomCoords = removeAtomsFarAway(atomCoords,maxParam,minParam)

%Array to store indices of rows to delete
deleteRowIndicesArray = zeros(length(atomCoords),1);

%Loop through all atoms
for eachAtom = 1 : length(atomCoords)
    
    %get (x,y,z) coordinates for the atom
    atomXYZ = cell2mat(atomCoords(eachAtom,9:11));
    
    %Check if the position of the atom is outside the box. If so then
    %remove the atom
    if atomXYZ(1) < minParam(1) || atomXYZ(1) > maxParam(1)
        deleteRowIndicesArray(eachAtom) = eachAtom;
    elseif atomXYZ(2) < minParam(2) || atomXYZ(2) > maxParam(2)
        deleteRowIndicesArray(eachAtom) = eachAtom;
    elseif atomXYZ(3) < minParam(3) || atomXYZ(3) > maxParam(3)
        deleteRowIndicesArray(eachAtom) = eachAtom;
    end
end

%find indices where element value = 0 corresponding to atoms we want to
%keep
keepRowIndex = ~deleteRowIndicesArray;

%transfer data for kept atoms to new cell array
newAtomCoords = atomCoords(keepRowIndex,:);

end

