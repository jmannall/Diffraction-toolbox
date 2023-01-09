clear all

n = 1;
wedgeIndex = 180:n:360;
minAngle = 0:n:180;
bendingAngle = n:n:360;

numObs = 20000;

geometry = RandomGeometryWedge(numObs);

w = geometry.wedgeIndex;
bA = geometry.bendingAngle;
mA = geometry.minAngle;

test = [w, bA, mA];
num = length(test);

%% Plots

close all

figure
h = histogram(w, length(wedgeIndex));
title('wedgeIndex')

figure
h = histogram(bA);
title('bendingAngle')


figure
h = histogram(mA);
title('minAngle')
