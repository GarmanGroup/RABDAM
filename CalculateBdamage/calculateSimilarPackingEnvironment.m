%% Function to group atoms into groups of similar packing density environments
function atomCoordsWithGroupings = calculateSimilarPackingEnvironment(atomCoords,binSize)

atomCoordsWithGroupings = atomCoords;

%Get the packing densities from the atom coordinate information
packingDensities = cell2mat(atomCoords(:,18));

%Calculate the minimum and maximum packing density values
maxPackingDensity = max(packingDensities);
minPackingDensity = min(packingDensities);

%calculate the minimum and maximum values for the loop
startValue = floor(minPackingDensity/binSize) * binSize;
endValue = floor(maxPackingDensity/binSize) * binSize;

%Create the array that the loop will go through 
arrayForLoop = startValue : binSize : endValue;

%Create a variable that will record the number of iterations. Each
%iteration number will correspond to an similar environment group
similarEnvironmentNumber = 0;

%Loop through each bin 
for binGroup = arrayForLoop
    %Increment the similar environment group number
    similarEnvironmentNumber = similarEnvironmentNumber + 1;
    
    %Calculate Environment variable for this group
    comparisonValue = floor(binGroup/binSize);
    
    %Calculate the binning condition for current similarity group
    indicesForGroup = find(floor(packingDensities./binSize) == comparisonValue);
    
    %Number of atoms in similar environment
    numberOfAtomsInEnvironment = length(indicesForGroup);
    
    for eachGroupMember = 1 : numberOfAtomsInEnvironment
        %Store the similar environment group number 
        atomCoordsWithGroupings{indicesForGroup(eachGroupMember),19} = similarEnvironmentNumber;
        %Store the explicit atomic number bin definition
        atomCoordsWithGroupings{indicesForGroup(eachGroupMember),20} = sprintf('%.0f <= x < %.0f',binGroup,binGroup + binSize);
        %Store the number of atoms in this similar packing density
        %environment
        atomCoordsWithGroupings{indicesForGroup(eachGroupMember),21} = numberOfAtomsInEnvironment;
    end
end

end

