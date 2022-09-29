%% Create result and geometry structs

close all

set(groot, 'defaultLineLineWidth', 1)


wedgeLength = 20;
radiusS = 1;
radiusR = 1;
zS = 10;
zR = 10;
fs = 96000;
wedge = 190:10:350;
bendingAngle = 190:10:350;
minAngle = 0:10:180;

[result, geometry] = SingleWedgeArray(wedgeLength, radiusS, radiusR, zS, zR, fs, wedge, bendingAngle, minAngle);

result = rmfield(result,'i');

numResults = length(result);
count = 0;

%% Process
input = [[geometry.wedge], [geometry.bendingAngle], [geometry.minAngle]];

[wedge, vw, iw] = unique(input(:,1));
[bendingAngle, vbA, ibA] = unique(input(:,2));
[minAngle, vmA, imA] = unique(input(:,3));

numWedges = length(wedge);
numBendingAngles = length(bendingAngle);
numMinAngles = length(minAngle);

%% Wedge

ibAimA = [ibA, imA];
[bAmA, vbAmA, ibAmA] = unique(ibAimA, 'rows');
numbAmA = length(bAmA);

control = struct('c1', bendingAngle, 'c2', minAngle, 'c1Name', 'Bending Angle', 'c2Name', 'Minimum Angle', 'cIndex', bAmA);
variable = struct('v', wedge, 'vName', 'Wedge');

% Save plots
ComparisonPlot(numbAmA, ibAmA, result, control, variable, true)
disp('Wedge figures complete');

%% Bending angle

iwimA = [iw, imA];
[wmA, vwmA, iwmA] = unique(iwimA, 'rows');
numwmA = length(wmA);

control = struct('c1', wedge, 'c2', minAngle, 'c1Name', 'Wedge', 'c2Name', 'Minimum Angle', 'cIndex', wmA);
variable = struct('v', bendingAngle, 'vName', 'Bending angle');

% Save plots
ComparisonPlot(numwmA, iwmA, result, control, variable, false)
disp('Bending angle figures complete');

%% Minimum angle

iwibA = [iw, ibA];
[wbA, vwbA, iwbA] = unique(iwibA, 'rows');
numwbA = length(wbA);

control = struct('c1', wedge, 'c2', bendingAngle, 'c1Name', 'Wedge', 'c2Name', 'Bending Angle', 'cIndex', wbA);
variable = struct('v', minAngle, 'vName', 'Minimum angle');

% Save plots
ComparisonPlot(numwbA, iwbA, result, control, variable, false)
disp('Minimum angle figures complete');

