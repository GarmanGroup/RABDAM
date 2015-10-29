%% Function to create input script for PDBCUR
function createInputFileForPDBCUR(fileName,inputs)

%Store all inputs in a cell array
allInputs = fieldnames(inputs);

%Open a text file for writing
fileID = fopen(fileName,'wt');

%Loop over all input lines
for inputLine = 1 : length(allInputs)
    %The first input is the filename so we can skip this line and only
    %write the lines that are not equal to the first line.
    if inputLine ~= 1
        %write line to input file
        fprintf(fileID,'%s\n',inputs.(allInputs{inputLine}));
    end
end

%close the file
fclose(fileID);

end

