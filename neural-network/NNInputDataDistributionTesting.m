clear all

n = 1;
wedgeIndex = 180:n:360;
minAngle = 0:n:180;
bendingAngle = n:n:360;

geometry = GeometryWedge(wedgeIndex, bendingAngle, minAngle, true, true);

w = geometry.wedgeIndex;
bA = geometry.bendingAngle;
mA = geometry.minAngle;

test = [w, bA, mA];
num = length(test);

%% Plots

close all

figure
h = histogram(w, length(wedgeIndex));
numBins = h.NumBins;
title('wedgeIndex')

pd = makedist('Triangular','A',180,'B',360,'C',360);
wDist = random(pd, num, 1);

figure
histogram(wDist, numBins)
title('wedgeIndex distribution')

figure
h = histogram(bA);
numBins = h.NumBins;
title('bendingAngle')

pd = makedist('Triangular','A',180,'B',180,'C',360);
bADist = random(pd, num, 1);

figure
histogram(bADist, numBins)
title('bendingAngle distribution')

figure
h = histogram(mA);
numBins = h.NumBins;
title('minAngle')

pd = makedist('Triangular','A',0,'B',0,'C',90);
mADist = random(pd, num, 1);

figure
histogram(mADist, numBins)
title('minAngle distribution')