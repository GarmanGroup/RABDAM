%% Function to calculate packing density of atoms 
function pdbAtomCoordsWithPackingDensity = calculatePackingDensityAlternative(pdbAtomCoords,allSurroundingAtomCoords,packingDensityThreshold)

%calculate the number of atoms in the PDB file
numberOfAtomsInPDB = length(pdbAtomCoords);

%calculate number of surrounding atoms
numberOfSurroundingAtoms = length(allSurroundingAtomCoords);

%Preallocate cell array to store packing density and other information
pdbAtomCoordsWithPackingDensity = cell(numberOfAtomsInPDB,21);

%put the atom coordinate information in the relevant elements of the new
%cell array
pdbAtomCoordsWithPackingDensity(:,1:15) = pdbAtomCoords;

%Get the x, y and z coordinates of the atoms in the original PDB file
xyzPDBAtom = str2double(pdbAtomCoords(:,9:11));

%Get the x, y and z coordinates of the surrounding atoms
xyzSurroundingAtom = cell2mat(allSurroundingAtomCoords(:,9:11));

%Function to loop through each atom and calculate the packing density for
%each one
packingDensityArray = calcPackingDensityLoop_mex(xyzPDBAtom,xyzSurroundingAtom,numberOfAtomsInPDB,numberOfSurroundingAtoms,packingDensityThreshold);

%Put the packing densities in the 18th column of the atom coordinates array
pdbAtomCoordsWithPackingDensity(:,18) = num2cell(packingDensityArray);

end