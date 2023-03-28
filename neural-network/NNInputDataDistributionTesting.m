clear all
%close all

% epsilon = 1e-6;
% wedgeIndex = [181, 360 - epsilon];
% 
% wI = linspace(wedgeIndex(1), wedgeIndex(2), 1e3);
% pdf = (1 / 8) * (wI - wedgeIndex(1)) .^ 4;
% cdf = (1 / 40) * (wI - wedgeIndex(1)) .^ 5;
% 
% b = 1 / cdf(end);
% cdf = (b / 40) * (wI - wedgeIndex(1)) .^ 5;
% 
% wI = nthroot((40 / b) * cdf, 5) + wedgeIndex(1);
% 
% x = RandomUniformDistribution([0 1], 20e3);
% test = nthroot(x, 5) * (wedgeIndex(2) - wedgeIndex(1)) + wedgeIndex(1);
% test2 = nthroot((40 / b) * x, 5) + wedgeIndex(1);
% 
% figure
% plot(wI, pdf)
% hold on
% plot(wI, cdf)
% 
% wI = RandomUniformDistribution([0 1], 20e3);
% wI = wedgeIndex(1) + (wedgeIndex(2) - wedgeIndex(1)) * 0.125 * (1 - wI .^ 4);
% 
% figure
% h = histogram(test);
% title('wedgeIndex')
% 
% figure
% h = histogram(test2);
% title('wedgeIndex2')

n = 1;
wedgeIndex = 180:n:360;
minAngle = 0:n:180;
bendingAngle = n:n:360;

numObs = 20000;

geometry = RandomGeometryWedge_Run2(numObs);

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

phii = atand(rS ./ abs(zA - zS));
%phii2 = atand(abs(zA - zR) ./ rR);


%% Plots

%close all

figure
h = histogram(w, length(wedgeIndex));
title('wedgeIndex')

figure
h = histogram(bA);
title('bendingAngle')

figure
h = histogram(mA);
title('minAngle')

figure
h = histogram(phii);
title('Phii')

figure
h = histogram(rS);
title('rS')

figure
h = histogram(rR);
title('rR')

L = sqrt((rS + rR) .^ 2 + abs(zS - zR) .^ 2);
[~,edges] = histcounts(log10(L));
figure
histogram(L,10.^edges)
set(gca, 'xscale','log')
title('L')

figure
h = histogram(L);
title('L')

figure
h = histogram(zS);
title('zS')

figure
h = histogram(zR);
title('zR')

