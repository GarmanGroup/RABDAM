%% Function to loop through the atoms in the PDB and determine their packing densities
function packingDensityArray = calcPackingDensityLoop(xyzPDBAtom,xyzSurroundingAtom,numberOfAtomsInPDB,numberOfSurroundingAtoms,packingDensityThreshold)

%Preallocating memory for the packingDensityArray
packingDensityArray = zeros(numberOfAtomsInPDB,1);

%Loop through each atom
for eachPDBAtom = 1 : numberOfAtomsInPDB

    %Set the atom packing density to zero. Here the packing density is
    %defined by the atomic contact number with contact being defined as any
    %atom within the packingDensityThreshold distance of it.
    atomPackingDensity = 0;

    %Loop through each of the surrounding atoms that can potentially
    %contribute to the packing density
    for eachSurroundingAtom = 1 : numberOfSurroundingAtoms

        %Calculate distance between the surrounding atom and the pdb atom
        distanceBetweenAtoms = sqrt((xyzPDBAtom(eachPDBAtom,1) - xyzSurroundingAtom(eachSurroundingAtom,1))^2 + (xyzPDBAtom(eachPDBAtom,2) - xyzSurroundingAtom(eachSurroundingAtom,2))^2 + (xyzPDBAtom(eachPDBAtom,3) - xyzSurroundingAtom(eachSurroundingAtom,3))^2);

        %Check if the distance is below the packing density threshold
        if distanceBetweenAtoms <= packingDensityThreshold
            %if so then increase the packing density counter by 1
            atomPackingDensity = atomPackingDensity + 1;
        end
    end
    %store the packing density of the atom in the cell array
    packingDensityArray(eachPDBAtom) = atomPackingDensity;
end


end
