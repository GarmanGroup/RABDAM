function newAtomCoords = removeAtomsFarAwayAlternative(atomCoords,maxParam,minParam)

%Get the atom coordinates from the atom coordinates array
atomXYZ = cell2mat(atomCoords(:,9:11));

%preallocate a logical array to hold values for position for each atom
conditionArray = false(length(atomCoords),3);

%Check conditions to find out whether the atom is inside or outside the
%box. If atom lies outside the x, y or z bounds then let the condition be
%for that coordinate be true
conditionArray(:,1) = atomXYZ(:,1) < minParam(1) | atomXYZ(:,1) > maxParam(1);
conditionArray(:,2) = atomXYZ(:,2) < minParam(2) | atomXYZ(:,2) > maxParam(2);
conditionArray(:,3) = atomXYZ(:,3) < minParam(3) | atomXYZ(:,3) > maxParam(3);

%Find the indices for the atoms that are inside the box. These correspond
%to the rows that contain all "false" values
keepRowIndex = ~any(conditionArray,2);

%Only store the information of the atoms that lie inside the box
newAtomCoords = atomCoords(keepRowIndex,:);

end

