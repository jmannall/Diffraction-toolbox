close all
controlparameters = struct('fs',48000, 'difforder', 1, 'nfft', 8192);
controlparameters.difforder = 3;

x = SinglePanel(0.5, 2, 2, controlparameters, true);

%% Test 2

wedgeIndex = 345;
minAngle = 20;
bendingAngle = 10:1:330;
geometry = GeometryWedge(wedgeIndex, bendingAngle, minAngle, false);
controlparameters = struct('fs',48000, 'difforder', 1, 'nfft', 8192);

[res, geom] = SingleWedgeArray(geometry, 10, 1, 1, 5, 5, controlparameters);