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


rS = geometry.rS;
rR = geometry.rR;
zS = geometry.zS;
zR = geometry.zR;
wedgeLength = geometry.wedgeLength;

zA = CalculateApex(rS, rR, zS, zR, wedgeLength, true);
zA = zA(:,3);

phii = 90 - atand(abs(zA - zS) ./ rS);
%phii2 = atand(abs(zA - zR) ./ rR);


%% Plots

close all

figure
h = histogram(w, length(wedgeIndex));
title('wedgeIndex')

figure
h = histogram(bA);
title('bendingAngle')

test = sum(mA < 15);
test2 = sum(mA >= 15);


figure
h = histogram(mA);
title('minAngle')

figure
h = histogram(phii);
title('Phii')

