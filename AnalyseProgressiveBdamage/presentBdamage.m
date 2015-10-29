%% Script to present Bdamage for all runs within a dataset
function presentBdamage(allBdamage,numberOfFiles,pathToBdamageFolder)

%Create puppet value for total point counter
totalPoints = 0;

%Find total datapoints across all sets
for file = 1 : numberOfFiles
    numberOfPoints = length(allBdamage{1,file});
    totalPoints = totalPoints + numberOfPoints;
end

%Create blank master x- and y-coordinate files
runNumber = zeros(totalPoints,1);
bdamage = zeros(totalPoints,1);

%Create a Point Counter and set value to zero
pointCounter = 0;

%Create a Max Bdamage Tracker and set value to zero
maxBdam = 0;

for file = 1 : numberOfFiles
    
    %Identify the number of atoms in the Bdamage file
    numberOfPoints = length(allBdamage{1,file});
    
    for point = 1 : numberOfPoints
        
        %Write run number to runNumber
        index = pointCounter + point;
        runNumber(index) = file;
        
        %Write Bdamage to bdamage
        bdamageValue = str2double(allBdamage{1,file}{point,17});
        bdamage(index) = bdamageValue;
        
        %Update the Max Bdamage Tracker if necessary
        if maxBdam < bdamageValue
            maxBdam = bdamageValue;
        end
    end
    
    %Update the point counter
    pointCounter = index;
end

%Create some buffer space at the top and RHS of the graph
xMax = numberOfFiles + 0.99;
yMax = (round(2*(maxBdam+0.5)))/2;

%Make a pretty scatter graph
fig = figure('name','Scatter Graph','Units', 'pixels','Position', [100 100 900 600]);
scatter(runNumber,bdamage,'x')
title('Bdamage Values for each Run','FontSize',20,'FontWeight','bold','FontName','AvantGarde')
xlabel('Dataset number','FontSize',18,'FontName','AvantGarde')
ylabel('Bdamage for each Atom','FontSize',18,'FontName','AvantGarde')

set(gca,                            ...
    'FontName'    , 'Helvetica'   , ...
    'FontSize'    , 14            , ...
    'Box'         , 'off'         , ...
    'TickDir'     , 'out'         , ...
    'TickLength'  , [.02 .02]     , ...
    'XLimMode'    , 'manual'      , ...
    'XLim'        , [0 xMax]      , ...
    'XTickMode'   , 'manual'      , ...
    'XTick'       , 1 : xMax      , ...
    'YLimMode'    , 'manual'      , ...
    'YLim'        , [0 yMax]      , ...
    'YMinorTick'  , 'on'          , ...
    'YGrid'       , 'on'          , ...
    'YMinorGrid'  , 'on'          , ...
    'XColor'      , [.3 .3 .3]    , ...
    'YColor'      , [.3 .3 .3]    , ...
    'LineWidth'   , 1);

%Print scatter graph to pdf
set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);
print(fig,'-dpdf', sprintf('%s\\%sBdamageProgression',pathToBdamageFolder,pathToBdamageFolder))

fprintf(sprintf('%sBdamageProgression.pdf created',pathToBdamageFolder))
fprintf('\n')

end