%% Script to present Bdamage for like atoms for two comparable datasets
function presentLikeAtoms(allLikeAtoms,pathToBdamageFolder)

%Find the number of Like Atoms
numberOfPoints = length(allLikeAtoms);
%numberOfPoints = 10505;

%Create a Max Bdamage Tracker and set value to zero
maxBdam = 0;

fprintf('AtomIDs of significant Bdamage ratio;\n')

%Create blank information cell
lineGraphData = zeros(numberOfPoints,2);

%Create blank Bdamage Ratio cell
bDamageRatio = zeros(numberOfPoints,1);
bDR2 = zeros(numberOfPoints,1);


for point = 1 : numberOfPoints
    
    %Write Bdamage to graph data file
    bdamageValue1 = str2double(allLikeAtoms{point,5});
    lineGraphData(point,1) = bdamageValue1;
    bdamageValue2 = str2double(allLikeAtoms{point,6});
    lineGraphData(point,2) = bdamageValue2;
    
    %calculate BdamageRatio
    ratio = lineGraphData(point,2)/lineGraphData(point,1);
    bDamageRatio(point,1) = (ratio + 1/ratio)/2;
    
    %calculate BdamageRatio2
    bD1 = lineGraphData(point,1);
    bD2 = lineGraphData(point,2);
    bDR2(point,1) = 2*(bD2-bD1)/(bD1+bD2);
    
    
    %return to user reference information of significant BdamageRatio
    if bDR2(point,1) > 0.25
        fprintf('%s %s %s %s\n',allLikeAtoms{point,1},allLikeAtoms{point,2},allLikeAtoms{point,3},allLikeAtoms{point,4})
    end
    
    %Update the max Bdamage Tracker if necessary
    if maxBdam < bdamageValue1
        maxBdam = bdamageValue1;
    end
    if maxBdam < bdamageValue2
        maxBdam = bdamageValue2;
    end
end

fprintf('\n')

%Create some buffer space at the top and RHS of the graph
xMax = 2.99;
yMax = (round(2*(maxBdam+0.5)))/2;

%Make a Kd plot
fig0 = figure('name','Distribution Plot','Units', 'pixels','Position', [100 100 900 600]);
for file = 1:2
pd = fitdist(lineGraphData(:,file),'Kernel','Kernel','epanechnikov');
x = 0:0.1:maxBdam;
y = pdf(pd,x);
plot(x,y)
hold on
end
hold off
title('\sl B\rm_{Damage} distribution plots')
xlabel('\sl B\rm_{Damage}')
ylabel('Frequency density')

set(gca,                            ...
    'FontName'    , 'Helvetica'   , ...
    'FontSize'    , 14            , ...
    'Box'         , 'off'         , ...
    'TickDir'     , 'out'         , ...
    'TickLength'  , [.02 .02]     , ...
    'XLimMode'    , 'manual'      , ...
    'XLim'        , [0 yMax]      , ...
    'XTickMode'   , 'manual'      , ...
    'XTick'       , 1 : yMax      , ...
    'YLimMode'    , 'auto'        , ...
    'YMinorTick'  , 'on'          , ...
    'YGrid'       , 'on'          , ...
    'YMinorGrid'  , 'on'          , ...
    'XColor'      , [.3 .3 .3]    , ...
    'YColor'      , [.3 .3 .3]    , ...
    'LineWidth'   , 1);

%Print distribution plot to pdf
set(fig0,'PaperOrientation','landscape');
set(fig0,'PaperUnits','normalized');
set(fig0,'PaperPosition', [0 0 1 1]);
print(fig0,'-dpdf', sprintf('%s\\%sBdamageDistribution',pathToBdamageFolder,pathToBdamageFolder))

fprintf(sprintf('%sBdamageDistribution.pdf created',pathToBdamageFolder))
fprintf('\n')

%make a boxplot
fig1 = figure('name','Scatter Graph','Units', 'pixels','Position', [100 100 900 600]);
boxplot(lineGraphData,'symbol','x')
title('\sl B\rm_{Damage} Values for each Dataset')
xlabel('Dataset number')
ylabel('\sl B\rm_{Damage} for each Atom')

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

%Print boxplot to pdf
set(fig1,'PaperOrientation','landscape');
set(fig1,'PaperUnits','normalized');
set(fig1,'PaperPosition', [0 0 1 1]);
print(fig1,'-dpdf', sprintf('%s\\%sLABdamageBPs',pathToBdamageFolder,pathToBdamageFolder))

fprintf(sprintf('%sLABdamageBPs.pdf created',pathToBdamageFolder))
fprintf('\n')

%make a stacked line graph
fig2 = figure('name','Line Graph','Units', 'pixels','Position', [100 100 900 600]);
plot(lineGraphData)
% hold on
% xa = [507,507];
% xb = [1039,1039];
% xc = [1611,1611];
% xd = [2099,2099];
% xe = [4247,4247];
% xf = [6349,6349];
% xg = [6842,6842];
% xh = [7350,7350];
% xi = [7889,7889];
% xj = [8375,8375];
% xk = [10506,10506];
% y = [0,yMax];
% plot(xa,y,'color','k','linewidth',1)
% hold on
% plot(xb,y,'color','k','linewidth',1)
% hold on
% plot(xc,y,'color','k','linewidth',1)
% hold on
% plot(xd,y,'color','k','linewidth',1)
% hold on
% plot(xe,y,'color','k','linewidth',1)
% hold on
% plot(xf,y,'color','k','linewidth',1)
% hold on
% plot(xg,y,'color','k','linewidth',1)
% hold on
% plot(xh,y,'color','k','linewidth',1)
% hold on
% plot(xi,y,'color','k','linewidth',1)
% hold on
% plot(xj,y,'color','k','linewidth',1)
% hold on
% plot(xk,y,'color','k','linewidth',1)
title('\sl B\rm_{Damage} Values for each Atom')
xlabel('Atom number')
ylabel('\sl B\rm_{Damage}')

set(gca,                            ...
    'FontName'    , 'Helvetica'   , ...
    'FontSize'    , 14            , ...
    'Box'         , 'off'         , ...
    'TickDir'     , 'out'         , ...
    'TickLength'  , [.02 .02]     , ...
    'XLimMode'    , 'auto'        , ...
    'XTickMode'   , 'auto'        , ...
    'YLimMode'    , 'manual'      , ...
    'YLim'        , [0 yMax]      , ...
    'YMinorTick'  , 'on'          , ...
    'YGrid'       , 'on'          , ...
    'YMinorGrid'  , 'on'          , ...
    'XColor'      , [.3 .3 .3]    , ...
    'YColor'      , [.3 .3 .3]    , ...
    'LineWidth'   , 1);

%Print stacked line graph to pdf
set(fig2,'PaperOrientation','landscape');
set(fig2,'PaperUnits','normalized');
set(fig2,'PaperPosition', [0 0 1 1]);
print(fig2,'-dpdf', sprintf('%s\\%sLABdamageSLG',pathToBdamageFolder,pathToBdamageFolder))

fprintf(sprintf('%sLABdamageSLG.pdf created',pathToBdamageFolder))
fprintf('\n')

%make a line graph for Bdamage ratio
fig3 = figure('name','Line Graph','Units', 'pixels','Position', [100 100 900 600]);
plot(bDamageRatio)
% hold on
% xa = [507,507];
% xb = [1039,1039];
% xc = [1611,1611];
% xd = [2099,2099];
% xe = [4247,4247];
% xf = [6349,6349];
% xg = [6842,6842];
% xh = [7350,7350];
% xi = [7889,7889];
% xj = [8375,8375];
% xk = [10506,10506];
% y = [0,1.04];
% plot(xa,y,'color','r','linewidth',1)
% hold on
% plot(xb,y,'color','r','linewidth',1)
% hold on
% plot(xc,y,'color','r','linewidth',1)
% hold on
% plot(xd,y,'color','r','linewidth',1)
% hold on
% plot(xe,y,'color','r','linewidth',1)
% hold on
% plot(xf,y,'color','r','linewidth',1)
% hold on
% plot(xg,y,'color','r','linewidth',1)
% hold on
% plot(xh,y,'color','r','linewidth',1)
% hold on
% plot(xi,y,'color','r','linewidth',1)
% hold on
% plot(xj,y,'color','r','linewidth',1)
% hold on
% plot(xk,y,'color','r','linewidth',1)
title('\sl B\rm_{DR} for each Atom')
xlabel('Atom number')
ylabel('\sl B\rm_{DR}')

set(gca,                            ...
    'FontName'    , 'Helvetica'   , ...
    'FontSize'    , 14            , ...
    'Box'         , 'off'         , ...
    'TickDir'     , 'out'         , ...
    'TickLength'  , [.02 .02]     , ...   
    'XLimMode'    , 'auto'        , ...
    'XTickMode'   , 'auto'        , ...
    'YLimMode'    , 'manual'      , ...
    'YLim'        , [1.01 1.04]   , ...
    'YMinorTick'  , 'on'          , ...
    'YGrid'       , 'on'          , ...
    'YMinorGrid'  , 'on'          , ...
    'XColor'      , [.3 .3 .3]    , ...
    'YColor'      , [.3 .3 .3]    , ...
    'LineWidth'   , 1);

%Print line graph to pdf
set(fig3,'PaperOrientation','landscape');
set(fig3,'PaperUnits','normalized');
set(fig3,'PaperPosition', [0 0 1 1]);
print(fig3,'-dpdf', sprintf('%s\\%sLABdamageRatio',pathToBdamageFolder,pathToBdamageFolder))

fprintf(sprintf('%sLABdamageRatio.pdf created',pathToBdamageFolder))
fprintf('\n')

%make a line graph for Bdamage ratio 2
fig4 = figure('name','Line Graph','Units', 'pixels','Position', [100 100 900 600]);
bar(bDR2)
% hold on
% xa = [507,507];
% xb = [1039,1039];
% xc = [1611,1611];
% xd = [2099,2099];
% xe = [4247,4247];
% xf = [6349,6349];
% xg = [6842,6842];
% xh = [7350,7350];
% xi = [7889,7889];
% xj = [8375,8375];
% xk = [10506,10505];
% y = [-0.2,0.3];
% plot(xa,y,'color','r','linewidth',1)
% hold on
% plot(xb,y,'color','r','linewidth',1)
% hold on
% plot(xc,y,'color','r','linewidth',1)
% hold on
% plot(xd,y,'color','r','linewidth',1)
% hold on
% plot(xe,y,'color','r','linewidth',1)
% hold on
% plot(xf,y,'color','r','linewidth',1)
% hold on
% plot(xg,y,'color','r','linewidth',1)
% hold on
% plot(xh,y,'color','r','linewidth',1)
% hold on
% plot(xi,y,'color','r','linewidth',1)
% hold on
% plot(xj,y,'color','r','linewidth',1)
% hold on
% plot(xk,y,'color','r','linewidth',1)
title('\sl B\rm_{DR} for each Atom')
xlabel('Atom number')
ylabel('\sl B\rm_{DR}')

set(gca,                            ...
    'FontName'    , 'Helvetica'   , ...
    'FontSize'    , 14            , ...
    'Box'         , 'off'         , ...
    'TickDir'     , 'out'         , ...
    'TickLength'  , [.02 .02]     , ...   
    'XLimMode'    , 'auto'        , ...
    'XTickMode'   , 'auto'        , ...
    'YLimMode'    , 'auto'        , ...
    'YMinorTick'  , 'on'          , ...
    'YGrid'       , 'on'          , ...
    'YMinorGrid'  , 'on'          , ...
    'XColor'      , [.3 .3 .3]    , ...
    'YColor'      , [.3 .3 .3]    , ...
    'LineWidth'   , 1);

%Print line graph to pdf
set(fig4,'PaperOrientation','landscape');
set(fig4,'PaperUnits','normalized');
set(fig4,'PaperPosition', [0 0 1 1]);
print(fig4,'-dpdf', sprintf('%s\\%sLABdamageRatio2',pathToBdamageFolder,pathToBdamageFolder))

fprintf(sprintf('%sLABdamageRatio2.pdf created',pathToBdamageFolder))
fprintf('\n')
    
end