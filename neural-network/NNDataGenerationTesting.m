close all
clear all

numObservations = 50e3;

geometry = RandomGeometryWedge(numObservations);

wedgeIndex = geometry.wedgeIndex;
bendingAngle = geometry.bendingAngle;
minAngle = geometry.minAngle;

figure
histogram(wedgeIndex)
title('wedgeIndex')

figure
histogram(bendingAngle)
title('bendingAngle')

figure
histogram(minAngle)
title('minAngle')
    