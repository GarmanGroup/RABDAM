%% Script to extract Bdamage values from only like atoms
%Input: Cell containing 2 parsed Bdamage files
function [allLikeAtoms] = parseLikeAtoms(allBdamage,stripHA,stripWaters)

%identify number of atoms in the first and second datasets 
numberOfPoints = length(allBdamage{1,1});
numberOfAtoms = length(allBdamage{1,2});

%create a puppet array into which information for like atoms can be
%shuttled
allLikeAtoms = cell(numberOfAtoms,6);

%set a row counter for output file
rowCounter = 1;

%if user has selected to ignore hetero atoms (implicitly including waters)
if stripHA
    %for each atom in file 1
    for point = 1 : numberOfPoints
        %ignore point if in a water molecule
        if ~strcmp(allBdamage{1,1}{point,1},'HETATM')
            %compare atom information to each atom in file 2
            for atom = 1 : numberOfAtoms
                %if the atom information is identical
                if strcmp(allBdamage{1,1}{point,3},allBdamage{1,2}{atom,3}) ...
                && strcmp(allBdamage{1,1}{point,5},allBdamage{1,2}{atom,5}) ...
                && strcmp(allBdamage{1,1}{point,6},allBdamage{1,2}{atom,6}) ...
                && strcmp(allBdamage{1,1}{point,7},allBdamage{1,2}{atom,7})
                    
                    %write the atom information and Bdamage Values to
                    %output array
                    allLikeAtoms{rowCounter,1} = allBdamage{1,1}{point,3};
                    allLikeAtoms{rowCounter,2} = allBdamage{1,1}{point,5};
                    allLikeAtoms{rowCounter,3} = allBdamage{1,1}{point,6};
                    allLikeAtoms{rowCounter,4} = allBdamage{1,1}{point,7};
                    allLikeAtoms{rowCounter,5} = allBdamage{1,1}{point,17};
                    allLikeAtoms{rowCounter,6} = allBdamage{1,2}{atom,17};
                    
                    %increment the row counter by 1
                    rowCounter = rowCounter + 1;
                end
            end
        end
    end
else
    %if user has selected to ignore waters
    if stripWaters
        %for each atom in file 1
        for point = 1 : numberOfPoints
            %ignore point if in a water molecule
            if ~strcmp(allBdamage{1,1}{point,5},'HOH')
%                 %compare atom information to each atom in file 2
%                 for atom = 1 : numberOfAtoms
%                     %if the atom information is identical
%                     if strcmp(allBdamage{1,1}{point,3},allBdamage{1,2}{atom,3}) ...
%                     && strcmp(allBdamage{1,1}{point,5},allBdamage{1,2}{atom,5}) ...
%                     && strcmp(allBdamage{1,1}{point,6},allBdamage{1,2}{atom,6}) ...
%                     && strcmp(allBdamage{1,1}{point,7},allBdamage{1,2}{atom,7}) ...
%                        
%                         %write the atom information and Bdamage Values to
%                         %output array
%                         allLikeAtoms{rowCounter,1} = allBdamage{1,1}{point,3};
%                         allLikeAtoms{rowCounter,2} = allBdamage{1,1}{point,5};
%                         allLikeAtoms{rowCounter,3} = allBdamage{1,1}{point,6};
%                         allLikeAtoms{rowCounter,4} = allBdamage{1,1}{point,7};
%                         allLikeAtoms{rowCounter,5} = allBdamage{1,1}{point,17};
%                         allLikeAtoms{rowCounter,6} = allBdamage{1,2}{atom,17};
%                         
%                         %increment the row counter by 1
%                         rowCounter = rowCounter + 1;
%                     end
%                 end
            end
        end
    %if user has not selected to ignore waters
    else
        %for each atom in file 1
        for point = 1 : numberOfPoints
            %compare atom information to each atom in file 2
            for atom = 1 : numberOfAtoms
                %if the atom information is identical
                if strcmp(allBdamage{1,1}{point,3},allBdamage{1,2}{atom,3}) ...
                && strcmp(allBdamage{1,1}{point,5},allBdamage{1,2}{atom,5}) ...
                && strcmp(allBdamage{1,1}{point,6},allBdamage{1,2}{atom,6}) ...
                && strcmp(allBdamage{1,1}{point,7},allBdamage{1,2}{atom,7}) ...
                    
                    %write the atom information and Bdamage Values to
                    %output array
                    allLikeAtoms{rowCounter,1} = allBdamage{1,1}{point,3};
                    allLikeAtoms{rowCounter,2} = allBdamage{1,1}{point,5};
                    allLikeAtoms{rowCounter,3} = allBdamage{1,1}{point,6};
                    allLikeAtoms{rowCounter,4} = allBdamage{1,1}{point,7};
                    allLikeAtoms{rowCounter,5} = allBdamage{1,1}{point,17};
                    allLikeAtoms{rowCounter,6} = allBdamage{1,2}{atom,17};
                   
                    %increment the row counter by 1
                    rowCounter = rowCounter + 1;
                end
            end
        end
    end
end

%remove empty rows from the resultant cell
allLikeAtoms(all(cellfun(@isempty,allLikeAtoms),2), : ) = [];

end