%% Function to discard the unit cells that do not contain any points in the box.

function remainingUnitCells = removeUnitCells(aTrans,bTrans,cTrans,atomCoords,maxParam,minParam)

%Preallocate memory for the unit cells that we will keep
remainingUnitCells = cell(1,27);
%% Set the coordinates of the opposite corners of the unit cell for each unit cell

%Match Unit Cell atom coords with their bounding unit cell coordinates
translateUnitCell = [0,0,0];
unitCell000 = cell(1,3);
unitCell000{1} = atomCoords.ac000; 
unitCell000{2} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);
translateUnitCell = [1,1,1];
unitCell000{3} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);

%Match Unit Cell atom coords with their bounding unit cell coordinates
translateUnitCell = [0,0,1];
unitCell001 = cell(1,3);
unitCell001{1} = atomCoords.ac001; 
unitCell001{2} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);
translateUnitCell = [1,1,2];
unitCell001{3} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);

%Match Unit Cell atom coords with their bounding unit cell coordinates
translateUnitCell = [0,0,-1];
unitCell002 = cell(1,3);
unitCell002{1} = atomCoords.ac002; 
unitCell002{2} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);
translateUnitCell = [1,1,0];
unitCell002{3} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);

%Match Unit Cell atom coords with their bounding unit cell coordinates
translateUnitCell = [0,1,0];
unitCell010 = cell(1,3);
unitCell010{1} = atomCoords.ac010;  
unitCell010{2} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);
translateUnitCell = [1,2,1];
unitCell010{3} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);

%Match Unit Cell atom coords with their bounding unit cell coordinates
translateUnitCell = [0,1,1];
unitCell011 = cell(1,3);
unitCell011{1} = atomCoords.ac011;  
unitCell011{2} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);
translateUnitCell = [1,2,2];
unitCell011{3} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);

%Match Unit Cell atom coords with their bounding unit cell coordinates
translateUnitCell = [0,1,-1];
unitCell012 = cell(1,3);
unitCell012{1} = atomCoords.ac012;  
unitCell012{2} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);
translateUnitCell = [1,2,0];
unitCell012{3} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);

%Match Unit Cell atom coords with their bounding unit cell coordinates
translateUnitCell = [0,-1,0];
unitCell020 = cell(1,3);
unitCell020{1} = atomCoords.ac020;  
unitCell020{2} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);
translateUnitCell = [1,0,1];
unitCell020{3} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);

%Match Unit Cell atom coords with their bounding unit cell coordinates
translateUnitCell = [0,-1,1];
unitCell021 = cell(1,3);
unitCell021{1} = atomCoords.ac021;  
unitCell021{2} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);
translateUnitCell = [1,0,2];
unitCell021{3} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);

%Match Unit Cell atom coords with their bounding unit cell coordinates
translateUnitCell = [0,-1,-1];
unitCell022 = cell(1,3);
unitCell022{1} = atomCoords.ac022;  
unitCell022{2} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);
translateUnitCell = [1,0,0];
unitCell022{3} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);

%Match Unit Cell atom coords with their bounding unit cell coordinates
translateUnitCell = [1,0,0];
unitCell100 = cell(1,3);
unitCell100{1} = atomCoords.ac100; 
unitCell100{2} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);
translateUnitCell = [2,1,1];
unitCell100{3} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);

%Match Unit Cell atom coords with their bounding unit cell coordinates
translateUnitCell = [1,0,1];
unitCell101 = cell(1,3);
unitCell101{1} = atomCoords.ac101; 
unitCell101{2} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);
translateUnitCell = [2,1,2];
unitCell101{3} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);

%Match Unit Cell atom coords with their bounding unit cell coordinates
translateUnitCell = [1,0,-1];
unitCell102 = cell(1,3);
unitCell102{1} = atomCoords.ac102; 
unitCell102{2} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);
translateUnitCell = [2,1,0];
unitCell102{3} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);

%Match Unit Cell atom coords with their bounding unit cell coordinates
translateUnitCell = [1,1,0];
unitCell110 = cell(1,3);
unitCell110{1} = atomCoords.ac110;  
unitCell110{2} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);
translateUnitCell = [2,2,1];
unitCell110{3} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);

%Match Unit Cell atom coords with their bounding unit cell coordinates
translateUnitCell = [1,1,1];
unitCell111 = cell(1,3);
unitCell111{1} = atomCoords.ac111;  
unitCell111{2} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);
translateUnitCell = [2,2,2];
unitCell111{3} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);

%Match Unit Cell atom coords with their bounding unit cell coordinates
translateUnitCell = [1,1,-1];
unitCell112 = cell(1,3);
unitCell112{1} = atomCoords.ac112;  
unitCell112{2} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);
translateUnitCell = [2,2,0];
unitCell112{3} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);

%Match Unit Cell atom coords with their bounding unit cell coordinates
translateUnitCell = [1,-1,0];
unitCell120 = cell(1,3);
unitCell120{1} = atomCoords.ac120;  
unitCell120{2} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);
translateUnitCell = [2,0,1];
unitCell120{3} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);

%Match Unit Cell atom coords with their bounding unit cell coordinates
translateUnitCell = [1,-1,1];
unitCell121 = cell(1,3);
unitCell121{1} = atomCoords.ac121;  
unitCell121{2} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);
translateUnitCell = [2,0,2];
unitCell121{3} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);

%Match Unit Cell atom coords with their bounding unit cell coordinates
translateUnitCell = [1,-1,-1];
unitCell122 = cell(1,3);
unitCell122{1} = atomCoords.ac122;  
unitCell122{2} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);
translateUnitCell = [2,0,0];
unitCell122{3} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);

%Match Unit Cell atom coords with their bounding unit cell coordinates
translateUnitCell = [-1,0,0];
unitCell200 = cell(1,3);
unitCell200{1} = atomCoords.ac200; 
unitCell200{2} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);
translateUnitCell = [0,1,1];
unitCell200{3} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);

%Match Unit Cell atom coords with their bounding unit cell coordinates
translateUnitCell = [-1,0,1];
unitCell201 = cell(1,3);
unitCell201{1} = atomCoords.ac201; 
unitCell201{2} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);
translateUnitCell = [0,1,2];
unitCell201{3} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);

%Match Unit Cell atom coords with their bounding unit cell coordinates
translateUnitCell = [-1,0,-1];
unitCell202 = cell(1,3);
unitCell202{1} = atomCoords.ac202; 
unitCell202{2} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);
translateUnitCell = [0,1,0];
unitCell202{3} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);

%Match Unit Cell atom coords with their bounding unit cell coordinates
translateUnitCell = [-1,1,0];
unitCell210 = cell(1,3);
unitCell210{1} = atomCoords.ac210;  
unitCell210{2} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);
translateUnitCell = [0,2,1];
unitCell210{3} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);

%Match Unit Cell atom coords with their bounding unit cell coordinates
translateUnitCell = [-1,1,1];
unitCell211 = cell(1,3);
unitCell211{1} = atomCoords.ac211;  
unitCell211{2} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);
translateUnitCell = [0,2,2];
unitCell211{3} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);

%Match Unit Cell atom coords with their bounding unit cell coordinates
translateUnitCell = [-1,1,-1];
unitCell212 = cell(1,3);
unitCell212{1} = atomCoords.ac212;  
unitCell212{2} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);
translateUnitCell = [0,2,0];
unitCell212{3} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);

%Match Unit Cell atom coords with their bounding unit cell coordinates
translateUnitCell = [-1,-1,0];
unitCell220 = cell(1,3);
unitCell220{1} = atomCoords.ac220;  
unitCell220{2} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);
translateUnitCell = [0,0,1];
unitCell220{3} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);

%Match Unit Cell atom coords with their bounding unit cell coordinates
translateUnitCell = [-1,-1,1];
unitCell221 = cell(1,3);
unitCell221{1} = atomCoords.ac221;  
unitCell221{2} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);
translateUnitCell = [0,0,2];
unitCell221{3} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);

%Match Unit Cell atom coords with their bounding unit cell coordinates
translateUnitCell = [-1,-1,-1];
unitCell222 = cell(1,3);
unitCell222{1} = atomCoords.ac222;  
unitCell222{2} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);
translateUnitCell = [0,0,0];
unitCell222{3} = aTrans*translateUnitCell(1) + bTrans*translateUnitCell(2) + cTrans*translateUnitCell(3);

%% Get corner points of asymmetric unit box

%Store the corner points of the box surrounding the atoms in the PDB
%coordinate file
corner1 = minParam;
corner2 = [minParam(1),minParam(2),maxParam(3)];
corner3 = [minParam(1),maxParam(2),minParam(3)];
corner4 = [minParam(1),maxParam(2),maxParam(3)];
corner5 = [maxParam(1),minParam(2),minParam(3)];
corner6 = [maxParam(1),minParam(2),maxParam(3)];
corner7 = [maxParam(1),maxParam(2),minParam(3)];
corner8 = maxParam;

%store the corner points in a 2d array
corners = [corner1;corner2;corner3;corner4;corner5;corner6;corner7;corner8];

%% Find the unit cells that we need to consider

%For each unit cell check whether any of the corner points of the box are
%contained within it. This is done for all 27 unit cells

for vertex = 1 : 8
   if (corners(vertex,1) >= unitCell000{2}(1) && corners(vertex,1) <= unitCell000{3}(1)...
           && corners(vertex,2) >= unitCell000{2}(2) && corners(vertex,2) <= unitCell000{3}(2)...
           && corners(vertex,3) >= unitCell000{2}(3) && corners(vertex,3) <= unitCell000{3}(3))
       remainingUnitCells{1} = unitCell000{1};
       break
   end
end

for vertex = 1 : 8
   if (corners(vertex,1) >= unitCell001{2}(1) && corners(vertex,1) <= unitCell001{3}(1)...
           && corners(vertex,2) >= unitCell001{2}(2) && corners(vertex,2) <= unitCell001{3}(2)...
           && corners(vertex,3) >= unitCell001{2}(3) && corners(vertex,3) <= unitCell001{3}(3))
       remainingUnitCells{2} = unitCell001{1};
       break
   end
end

for vertex = 1 : 8
   if (corners(vertex,1) >= unitCell002{2}(1) && corners(vertex,1) <= unitCell002{3}(1)...
           && corners(vertex,2) >= unitCell002{2}(2) && corners(vertex,2) <= unitCell002{3}(2)...
           && corners(vertex,3) >= unitCell002{2}(3) && corners(vertex,3) <= unitCell002{3}(3))
       remainingUnitCells{3} = unitCell002{1};
       break
   end
end

for vertex = 1 : 8
   if (corners(vertex,1) >= unitCell010{2}(1) && corners(vertex,1) <= unitCell010{3}(1)...
           && corners(vertex,2) >= unitCell010{2}(2) && corners(vertex,2) <= unitCell010{3}(2)...
           && corners(vertex,3) >= unitCell010{2}(3) && corners(vertex,3) <= unitCell010{3}(3))
       remainingUnitCells{4} = unitCell010{1};
       break
   end
end

for vertex = 1 : 8
   if (corners(vertex,1) >= unitCell011{2}(1) && corners(vertex,1) <= unitCell011{3}(1)...
           && corners(vertex,2) >= unitCell011{2}(2) && corners(vertex,2) <= unitCell011{3}(2)...
           && corners(vertex,3) >= unitCell011{2}(3) && corners(vertex,3) <= unitCell011{3}(3))
       remainingUnitCells{5} = unitCell011{1};
       break
   end
end

for vertex = 1 : 8
   if (corners(vertex,1) >= unitCell012{2}(1) && corners(vertex,1) <= unitCell012{3}(1)...
           && corners(vertex,2) >= unitCell012{2}(2) && corners(vertex,2) <= unitCell012{3}(2)...
           && corners(vertex,3) >= unitCell012{2}(3) && corners(vertex,3) <= unitCell012{3}(3))
       remainingUnitCells{6} = unitCell012{1};
       break
   end
end

for vertex = 1 : 8
   if (corners(vertex,1) >= unitCell020{2}(1) && corners(vertex,1) <= unitCell020{3}(1)...
           && corners(vertex,2) >= unitCell020{2}(2) && corners(vertex,2) <= unitCell020{3}(2)...
           && corners(vertex,3) >= unitCell020{2}(3) && corners(vertex,3) <= unitCell020{3}(3))
       remainingUnitCells{7} = unitCell020{1};
       break
   end
end

for vertex = 1 : 8
   if (corners(vertex,1) >= unitCell021{2}(1) && corners(vertex,1) <= unitCell021{3}(1)...
           && corners(vertex,2) >= unitCell021{2}(2) && corners(vertex,2) <= unitCell021{3}(2)...
           && corners(vertex,3) >= unitCell021{2}(3) && corners(vertex,3) <= unitCell021{3}(3))
       remainingUnitCells{8} = unitCell021{1};
       break
   end
end

for vertex = 1 : 8
   if (corners(vertex,1) >= unitCell022{2}(1) && corners(vertex,1) <= unitCell022{3}(1)...
           && corners(vertex,2) >= unitCell022{2}(2) && corners(vertex,2) <= unitCell022{3}(2)...
           && corners(vertex,3) >= unitCell022{2}(3) && corners(vertex,3) <= unitCell022{3}(3))
       remainingUnitCells{9} = unitCell022{1};
       break
   end
end

for vertex = 1 : 8
   if (corners(vertex,1) >= unitCell100{2}(1) && corners(vertex,1) <= unitCell100{3}(1)...
           && corners(vertex,2) >= unitCell100{2}(2) && corners(vertex,2) <= unitCell100{3}(2)...
           && corners(vertex,3) >= unitCell100{2}(3) && corners(vertex,3) <= unitCell100{3}(3))
       remainingUnitCells{10} = unitCell100{1};
       break
   end
end

for vertex = 1 : 8
   if (corners(vertex,1) >= unitCell101{2}(1) && corners(vertex,1) <= unitCell101{3}(1)...
           && corners(vertex,2) >= unitCell101{2}(2) && corners(vertex,2) <= unitCell101{3}(2)...
           && corners(vertex,3) >= unitCell101{2}(3) && corners(vertex,3) <= unitCell101{3}(3))
       remainingUnitCells{11} = unitCell101{1};
       break
   end
end

for vertex = 1 : 8
   if (corners(vertex,1) >= unitCell102{2}(1) && corners(vertex,1) <= unitCell102{3}(1)...
           && corners(vertex,2) >= unitCell102{2}(2) && corners(vertex,2) <= unitCell102{3}(2)...
           && corners(vertex,3) >= unitCell102{2}(3) && corners(vertex,3) <= unitCell102{3}(3))
       remainingUnitCells{12} = unitCell102{1};
       break
   end
end

for vertex = 1 : 8
   if (corners(vertex,1) >= unitCell110{2}(1) && corners(vertex,1) <= unitCell110{3}(1)...
           && corners(vertex,2) >= unitCell110{2}(2) && corners(vertex,2) <= unitCell110{3}(2)...
           && corners(vertex,3) >= unitCell110{2}(3) && corners(vertex,3) <= unitCell110{3}(3))
       remainingUnitCells{13} = unitCell110{1};
       break
   end
end

for vertex = 1 : 8
   if (corners(vertex,1) >= unitCell111{2}(1) && corners(vertex,1) <= unitCell111{3}(1)...
           && corners(vertex,2) >= unitCell111{2}(2) && corners(vertex,2) <= unitCell111{3}(2)...
           && corners(vertex,3) >= unitCell111{2}(3) && corners(vertex,3) <= unitCell111{3}(3))
       remainingUnitCells{14} = unitCell111{1};
       break
   end
end

for vertex = 1 : 8
   if (corners(vertex,1) >= unitCell112{2}(1) && corners(vertex,1) <= unitCell112{3}(1)...
           && corners(vertex,2) >= unitCell112{2}(2) && corners(vertex,2) <= unitCell112{3}(2)...
           && corners(vertex,3) >= unitCell112{2}(3) && corners(vertex,3) <= unitCell112{3}(3))
       remainingUnitCells{15} = unitCell112{1};
       break
   end
end

for vertex = 1 : 8
   if (corners(vertex,1) >= unitCell120{2}(1) && corners(vertex,1) <= unitCell120{3}(1)...
           && corners(vertex,2) >= unitCell120{2}(2) && corners(vertex,2) <= unitCell120{3}(2)...
           && corners(vertex,3) >= unitCell120{2}(3) && corners(vertex,3) <= unitCell120{3}(3))
       remainingUnitCells{16} = unitCell120{1};
       break
   end
end

for vertex = 1 : 8
   if (corners(vertex,1) >= unitCell121{2}(1) && corners(vertex,1) <= unitCell121{3}(1)...
           && corners(vertex,2) >= unitCell121{2}(2) && corners(vertex,2) <= unitCell121{3}(2)...
           && corners(vertex,3) >= unitCell121{2}(3) && corners(vertex,3) <= unitCell121{3}(3))
       remainingUnitCells{17} = unitCell121{1};
       break
   end
end

for vertex = 1 : 8
   if (corners(vertex,1) >= unitCell122{2}(1) && corners(vertex,1) <= unitCell122{3}(1)...
           && corners(vertex,2) >= unitCell122{2}(2) && corners(vertex,2) <= unitCell122{3}(2)...
           && corners(vertex,3) >= unitCell122{2}(3) && corners(vertex,3) <= unitCell122{3}(3))
       remainingUnitCells{18} = unitCell122{1};
       break
   end
end

for vertex = 1 : 8
   if (corners(vertex,1) >= unitCell200{2}(1) && corners(vertex,1) <= unitCell200{3}(1)...
           && corners(vertex,2) >= unitCell200{2}(2) && corners(vertex,2) <= unitCell200{3}(2)...
           && corners(vertex,3) >= unitCell200{2}(3) && corners(vertex,3) <= unitCell200{3}(3))
       remainingUnitCells{19} = unitCell200{1};
       break
   end
end

for vertex = 1 : 8
   if (corners(vertex,1) >= unitCell201{2}(1) && corners(vertex,1) <= unitCell201{3}(1)...
           && corners(vertex,2) >= unitCell201{2}(2) && corners(vertex,2) <= unitCell201{3}(2)...
           && corners(vertex,3) >= unitCell201{2}(3) && corners(vertex,3) <= unitCell201{3}(3))
       remainingUnitCells{20} = unitCell201{1};
       break
   end
end

for vertex = 1 : 8
   if (corners(vertex,1) >= unitCell202{2}(1) && corners(vertex,1) <= unitCell202{3}(1)...
           && corners(vertex,2) >= unitCell202{2}(2) && corners(vertex,2) <= unitCell202{3}(2)...
           && corners(vertex,3) >= unitCell202{2}(3) && corners(vertex,3) <= unitCell202{3}(3))
       remainingUnitCells{21} = unitCell202{1};
       break
   end
end

for vertex = 1 : 8
   if (corners(vertex,1) >= unitCell210{2}(1) && corners(vertex,1) <= unitCell210{3}(1)...
           && corners(vertex,2) >= unitCell210{2}(2) && corners(vertex,2) <= unitCell210{3}(2)...
           && corners(vertex,3) >= unitCell210{2}(3) && corners(vertex,3) <= unitCell210{3}(3))
       remainingUnitCells{22} = unitCell210{1};
       break
   end
end

for vertex = 1 : 8
   if (corners(vertex,1) >= unitCell211{2}(1) && corners(vertex,1) <= unitCell211{3}(1)...
           && corners(vertex,2) >= unitCell211{2}(2) && corners(vertex,2) <= unitCell211{3}(2)...
           && corners(vertex,3) >= unitCell211{2}(3) && corners(vertex,3) <= unitCell211{3}(3))
       remainingUnitCells{23} = unitCell211{1};
       break
   end
end

for vertex = 1 : 8
   if (corners(vertex,1) >= unitCell212{2}(1) && corners(vertex,1) <= unitCell212{3}(1)...
           && corners(vertex,2) >= unitCell212{2}(2) && corners(vertex,2) <= unitCell212{3}(2)...
           && corners(vertex,3) >= unitCell212{2}(3) && corners(vertex,3) <= unitCell212{3}(3))
       remainingUnitCells{24} = unitCell212{1};
       break
   end
end

for vertex = 1 : 8
   if (corners(vertex,1) >= unitCell220{2}(1) && corners(vertex,1) <= unitCell220{3}(1)...
           && corners(vertex,2) >= unitCell220{2}(2) && corners(vertex,2) <= unitCell220{3}(2)...
           && corners(vertex,3) >= unitCell220{2}(3) && corners(vertex,3) <= unitCell220{3}(3))
       remainingUnitCells{25} = unitCell220{1};
       break
   end
end

for vertex = 1 : 8
   if (corners(vertex,1) >= unitCell221{2}(1) && corners(vertex,1) <= unitCell221{3}(1)...
           && corners(vertex,2) >= unitCell221{2}(2) && corners(vertex,2) <= unitCell221{3}(2)...
           && corners(vertex,3) >= unitCell221{2}(3) && corners(vertex,3) <= unitCell221{3}(3))
       remainingUnitCells{26} = unitCell221{1};
       break
   end
end

for vertex = 1 : 8
   if (corners(vertex,1) >= unitCell222{2}(1) && corners(vertex,1) <= unitCell222{3}(1)...
           && corners(vertex,2) >= unitCell222{2}(2) && corners(vertex,2) <= unitCell222{3}(2)...
           && corners(vertex,3) >= unitCell222{2}(3) && corners(vertex,3) <= unitCell222{3}(3))
       remainingUnitCells{27} = unitCell222{1};
       break
   end
end

%Remove empty elements
remainingUnitCells = remainingUnitCells(~cellfun('isempty', remainingUnitCells));

end