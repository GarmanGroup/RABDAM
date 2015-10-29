%% Function to calculate BDamage for each atom
function allAtomCoords = calcBDamage(atomCoords)

allAtomCoords = atomCoords;

%Get the similar packing density environment group numbers
groupNumbers = cell2mat(atomCoords(:,19));

%Get the unique group number
maxGroupNumbers = max(groupNumbers);

for group = 1 : maxGroupNumbers
    %Find the indices for the rows that correspond to atoms in the cuurent
    %group
    rowIndicesForGroup = groupNumbers == group;
    
    %Get find the average bfactor for the atoms in the group
    averageBFactor = mean(str2double(atomCoords(rowIndicesForGroup,13)));
    
    %Convert from logical values for row indices to actual double values
    rowIndicesForGroup = find(rowIndicesForGroup);
    
    for eachAtom = 1 : length(rowIndicesForGroup)
        %Insert Average Bfactor for group into the array
        allAtomCoords{rowIndicesForGroup(eachAtom),16} = averageBFactor;
        
        %Get B factor for atom
        atomBFactor = str2double(atomCoords{rowIndicesForGroup(eachAtom),13});
        %Calculate BDamage and insert it into the array
        allAtomCoords{rowIndicesForGroup(eachAtom),17} = atomBFactor/averageBFactor;
    end
end

end

