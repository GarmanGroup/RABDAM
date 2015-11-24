%% Script to calculate Bprogression
function [allBdamage] = calculateBprogression(allBdamage,numberOfFiles,numberOfAtoms,pathToBdamageFolder)

%Calculate Bprogression for each file
for file = 1 : numberOfFiles
    
    %Identify the number of rows in the file
    param = size(allBdamage{1,file});
    numberOfRows = param(1);
    
    %Create puppet cell array into which Bprogression values will be
    %written
    bprogCol = cell(numberOfRows,1);
    
    %Calculate Bprogression for each atom
    for atom = 1 : numberOfAtoms
        
        %Calculate the value of Bprogression by taking the Bfactor for the
        %atom and dividing it by the average Bfactor of the class in which
        %the atom exists in the reference (first) data set
        bprog = ((str2double(allBdamage{1,file}{atom,13}))/(str2double(allBdamage{1,1}{atom,16})));
        
        %Write Bprogression to the coresponding row in the puppet cell
        %array
        bprogCol{atom} = bprog;
    end
    
    %Write the puppet cell array to be the final column of the
    %file's original extracted data
    allBdamage{1,file} = [allBdamage{1,file} bprogCol];  
end

%Write the Bdamage to output folder in the input directory
fprintf(sprintf('writing the calculated Bprogression values to output file %s\\ Bprogression.txt',pathToBdamageFolder))
fprintf('\n')

end