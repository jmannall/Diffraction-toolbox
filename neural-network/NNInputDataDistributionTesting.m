clear all
close all

fs = 96e3;
nfft = 16384;
c = 344;
controlparameters = struct('fs', fs, 'nfft', nfft, 'difforder', 1, 'c', c, 'saveFiles', 2, 'noDirect', true);
    
out = SingleWedge(10, 359.99, 0.1, 180.5, 1, 1, 5, 5, controlparameters, true);


w = [180 360];
mA = [0 (w(2) - 180) / 2];
bA = [w(1) - 2 * mA(1) w(2) - 2 * mA(1)];

corners = [bA(1) mA(1) w(1)
    bA(1) mA(1) w(2)
    bA(1) mA(2) w(2)
    bA(2) mA(1) w(2)];

planes = [1 4 2
    1 2 3
    2 4 3
    1 3 4];

figure
patch('Faces', planes, 'Vertices', corners, 'facealpha', 0.1)
view(3)

mAMax = 90;
bAMax = 360;
bAMin = 180;

wMax = 360;
wMin = 180;

base = mAMax * (bAMax - bAMin) / 2;
height = wMax - wMin;

volume = (1 / 3) * base * height;

scale = 180^3 / 12;
w = 180:1:360;
pdf = (1 / (4 * volume) * (w - wMin) .^ 2);
pdf2 = (1 / (60 * 180 ^ 2) * (w - 180) .^ 2);
cdf = (1/ (12 * volume) * (w - wMin) .^ 3);
cdf2 = (1/ (180 ^ 3) * (w - wMin) .^ 3);


%volume = (1 / 3) * (90 * 180) / 2 * 180 = (90 * 180 ^ 2) / 6 = 15 * 180 ^ 2

wMin = 181;
%volume = (1 / 3) * (90 * 180) / 2 * 179 = (90 * 180 * 179) / 6 = 15 * 180 * 179
%pdf = (1 / (60 * 180 * 179) * (w - 181) .^ 2);
%cdf = (1/ (180 * 180 * 179) * (w - 181) .^ 3);


%%
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

%%

n = 1;
wedgeIndex = 180:n:360;
minAngle = 0:n:180;
bendingAngle = n:n:360;

numObs = 20000;
numObs = 20e3;

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

