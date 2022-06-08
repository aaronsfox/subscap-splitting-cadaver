%% This script creates 3d visuals of relevant data to demonstrate the main
%  effects observed in the results.

%% Set-up

%Create lists of independent variables
tendonLoads = [{'0N'}, {'20N'}, {'40N'}];

%Set colours for independent variables
viridisCol = viridis;
loadCols = [viridisCol(end,:); viridisCol(round(length(viridis)/2),:); viridisCol(1,:)];

%% Load data

%Get the data from the .mat file
load('..\\Results\\Processed\\alignedGenericSurfaces.mat');

%Set humerus coordinate system origin for plotting
humerusOrigin = [humerusCS.Xc(1), humerusCS.Xc(2), humerusCS.Xc(3)];

%Load the group data
data = readtable('..\\Results\\Processed\\groupData.csv');

%% Create visuals

%Overall effect of load on SS4 line of action in scapular plan

%Extract the relevant data
visData = data(strcmp(data.region, 'ss4') & ...
    strcmp(data.plane, 'SP'), :);

%Split across loads and calculate mean and SD
for loadVar = 1:length(tendonLoads)
    
    %Get the table
    calcData = visData(strcmp(visData.load, tendonLoads{loadVar}), :);
    
    %Calculate mean/SD and store in array
    loaVals(loadVar,1) = mean(calcData(:,'lineOfAction').lineOfAction);
    loaVals(loadVar,2) = std(calcData(:,'lineOfAction').lineOfAction);
    
    %Calculate +/- values from mean and SD
    loaVals(loadVar,3) = loaVals(loadVar,1) - loaVals(loadVar,2);
    loaVals(loadVar,4) = loaVals(loadVar,1) + loaVals(loadVar,2);    
    
end

%Create base surfaces
cFigure; hold on

%Bones
gpatch(scapulaF,scapulaV,'kw','none',0.3);
gpatch(humerusF,humerusV,'kw','none',0.3);
axisGeom;
view(90,0); %scapula view; needs extra 90 rotation
camlight headlight;

%Set axis limits
ax = gca();
ax.YLim = [-125, ax.YLim(2)];
ax.Clipping = true;
axis off;

%Coordinate system (Y/Z of humerus system)
plotCS = [{'Yc'},{'Zc'}];
for cc = 1:length(plotCS)
    %Get points and vectors
    x = humerusCS.(plotCS{cc})(1); y = humerusCS.(plotCS{cc})(2); z = humerusCS.(plotCS{cc})(3);
    u = humerusCS.(plotCS{cc})(4); v = humerusCS.(plotCS{cc})(5); w = humerusCS.(plotCS{cc})(6);
    %Normalise lengths
    Ln = sqrt(u.^2 + v.^2 + w.^2);
    u = u./Ln; v = v./Ln; w = w./Ln;  %normalize vectors
    MaxLen = 1e-1*1000; %max length preference
    %Set vector length as max length
    u = u*MaxLen;
    v = v*MaxLen;  
    w = w*MaxLen;
    %Plot axes
    quiver3(x,y,z,u,v,w,'k','LineWidth',2,'ShowArrowHead','off')
end
clear cc

%Add mean and SD values for region with each load

%Set vector length
vecLength = 1e-1*1000;

%Loop through loads
for loadVar = 1:length(tendonLoads)
    
    %Get the mean value
    loa = loaVals(loadVar,1);
    
    %Mean line
    if loa < 180
        %Calculate coordinates
        vecX = 0;
        vecY = sind((180-loa)) * vecLength * -1;
        vecZ = cosd((180-loa)) * vecLength * -1;
        %Plot line of action    
        vLine = plot3([humerusOrigin(1),vecX],...
            [humerusOrigin(2),vecY],...
            [humerusOrigin(3),vecZ],...
            'Color',loadCols(loadVar,:),'LineWidth',2);
    elseif loa > 180
        %Calculate coordinates
        vecX = 0;
        vecY = sind((loa-180)) * vecLength;
        vecZ = cosd((loa-180)) * vecLength * -1;
        %Plot line of action    
        vLine = plot3([humerusOrigin(1),vecX],...
            [humerusOrigin(2),vecY],...
            [humerusOrigin(3),vecZ],...
            'Color',loadCols(loadVar,:),'LineWidth',2);
    end
    
end







