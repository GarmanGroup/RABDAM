%% Function to calculate packing density of atoms 
function pdbAtomCoordsWithPackingDensity = calculatePackingDensity(pdbAtomCoords,allSurroundingAtomCoords,packingDensityThreshold)

%calculate the number of atoms in the PDB file
numberOfAtomsInPDB = length(pdbAtomCoords);

%calculate number of surrounding atoms
numberOfSurroundingAtoms = length(allSurroundingAtomCoords);

%Preallocate cell array to store packing density and other information
pdbAtomCoordsWithPackingDensity = cell(numberOfAtomsInPDB,21);

%put the atom coordinate information in the relevant elements of the new
%cell array
pdbAtomCoordsWithPackingDensity(:,1:15) = pdbAtomCoords;

%Loop through each atom
for eachPDBAtom = 1 : numberOfAtomsInPDB 
    
    %Get the atom coordinates
    xyzPDBAtom = str2double(pdbAtomCoords(eachPDBAtom,9:11));
    
    %Set the atom packing density to zero. Here the packing density is
    %defined by the atomic contact number with contact being defined as any
    %atom within the packingDensityThreshold distance of it.
    atomPackingDensity = 0;
    
    %Loop through each of the surrounding atoms that can potentially
    %contribute to the packing density
    for eachSurroundingAtom = 1 : numberOfSurroundingAtoms
        %Get coordinates of surrounding atom
        xyzSurroundingAtom = cell2mat(allSurroundingAtomCoords(eachSurroundingAtom,9:11));
        
        %Calculate distance between the surrounding atom and the pdb atom
        distanceBetweenAtoms = pdist2(xyzPDBAtom,xyzSurroundingAtom);
        
        %Check if the distance is below the packing density threshold
        if distanceBetweenAtoms <= packingDensityThreshold
            %if so then increase the packing density counter by 1
            atomPackingDensity = atomPackingDensity + 1;
        end
    end
    %store the packing density of the atom in the cell array
    pdbAtomCoordsWithPackingDensity{eachPDBAtom,18} = atomPackingDensity;   
end

end