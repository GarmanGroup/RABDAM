%% Function to create a pdb file given the start, atomic coordinates and end of the file
function createPDBFile(atomicCoordinates, pdbPreamble, pdbEOF, fileName, fileNameSuffix)

%Split the file name because that contains ".pdb" at the end and we don't
%want that here
fileNameSplit = strsplit(fileName,{'.'});

%Take the first part of the split file name because that contains the part
%of the filename that we want
fileNamePrefix = fileNameSplit{1};

%open a file for writing the pdb entry and give it the specified filename.
fileID = fopen(sprintf('%s%s.pdb',fileNamePrefix,fileNameSuffix),'w');

%write the pdb preamble into the file
fprintf(fileID,'%s',pdbPreamble);

%loop over every atom
for eachAtomRecord = 1 : length(atomicCoordinates)
    %write the atom record into the pdb file
    fprintf(fileID,'%6s%5s %s%s%s %s%s%s%11.3f%8.3f%8.3f%s%s%12s%s\n',atomicCoordinates{eachAtomRecord,:});
end

%write the end of the pdb file
fprintf(fileID,'%s',pdbEOF);

%Close the file
fclose(fileID);

end

