close all

%% Data

fs = 96e3;
nfft = 4096;
numEdges = 2;

thetaS = 10;
thetaR = 260;
wedgeIndex = [270; 270];
W = 2;
radiusS = 0.5;
radiusR = 0.5;
height = 10;
[zS, zR] = deal(height / 2);

thetaS = 45;
thetaR = 205;
W = 0.89;
radiusS = 0.8;
radiusR = 0.23;

controlparameters = struct('fs', fs, 'nfft', nfft, 'difforder', numEdges);

%% Generate geometry

% data = struct('rS', radiusS, 'rR', radiusR, 'W', W, 'L', radiusS + sum(W) + radiusR, 'thetaS', thetaS, 'thetaR', thetaR, 'wedgeIndex', wedgeIndex);
% [source, receiver, Q, apex, corners, planeCorners, planeRigid, valid] = CreateNthOrderPathData(wedgeIndex, thetaS, thetaR, radiusS, radiusR, W, height);

[source, receiver, Q, apex, corners, planeCorners, planeRigid, data] = GenerateNthOrderPath(numEdges, height);


% Generate CAD file from source receiver and top down veiw
% Calculate geometry parameters (or add as output to above function) for
% 1st order diffraction sections.

% Expand to include z variation. Requires calculating all the apex points.

createPlot = false;
%% 2nd order BTM

[ir, tfmag, tvec, fvec, tfcomplex] = SingleBTM(source, receiver, corners, planeCorners, planeRigid, controlparameters, createPlot);

%% BTM Daisy Chains

withCorrection = true;
[apexTfmag, fvec, apexTfcomplex] = SingleBTMApexDaisyChain(source, receiver, apex, corners, planeCorners, controlparameters, data, withCorrection, createPlot);
[planeApexTfmag, fvec, planeApexTfcomplex] = SingleBTMPlaneDaisyChain(source, receiver, apex, corners, planeCorners, controlparameters, data, false, withCorrection, createPlot);
[planeTfmag, fvec, planeTfcomplex] = SingleBTMPlaneDaisyChain(source, receiver, apex, corners, planeCorners, controlparameters, data, true, withCorrection, createPlot);

withCorrection = false;
[apexTfmag2, fvec, apexTfcomplex2] = SingleBTMApexDaisyChain(source, receiver, apex, corners, planeCorners, controlparameters, data, withCorrection, createPlot);
[planeApexTfmag2, fvec, planeApexTfcomplex2] = SingleBTMPlaneDaisyChain(source, receiver, apex, corners, planeCorners, controlparameters, data, false, withCorrection, createPlot);
[planeTfmag2, fvec, planeTfcomplex2] = SingleBTMPlaneDaisyChain(source, receiver, apex, corners, planeCorners, controlparameters, data, true, withCorrection, createPlot);

%% UTD

always = true;
phi = asin(1);
c = 344;

withCorrection = false;
[utdTfmag, fvec, utdTfcomplex] = SingleUTDApexDaisyChain(data, phi, always, c, controlparameters, withCorrection);

% c = 344;
% wedgeIndex = 340;
% v = pi / deg2rad(wedgeIndex);
% thetaS = 10;
% radiusS = 1;
% thetaR = 300;
% radiusR = 1;
% 
% 
% 
% [tfmagUTD, fvec, tfcomplexUTD] = SingleUTDWedge(thetaS, thetaR, radiusS, radiusR, wedgeIndex, phi, always, controlparameters, c);
% 
% [ir, tfmag] = SingleWedge(20, wedgeIndex, thetaS, thetaR, radiusS, radiusR, 10, 10, controlparameters, false);
% 

%% Figures
% figure
% semilogx(fvec, tfmag.diff1)
% hold on
% semilogx(fvec, tfmagUTD)
% legend('BTM', 'UTD')


% figure('Position',[100, 50, 800, 400])
% semilogx(fvec, apexTfmag(:,end))
% hold on
% semilogx(fvec, planeTfmag(:,end))
% semilogx(fvec, planeApexTfmag(:,end))
% semilogx(fvec, utdTfmag(:,end))
% semilogx(fvec, tfmag.diff2)
% title('With correction')
% legend('BTM Apex Daisy Chain', 'BTM Plane Apex Daisy Chain', 'BTM Plane Daisy Chain', 'UTD', 'True BTM', 'Location', 'southwest')
% xlim([20 20000])
% ylim([-70 0])

figure('Position',[50 100 1820 800])
t = tiledlayout(1,5);
nexttile([1 2])
plot3(source(:,1), source(:,2), source(:,3), 'o')
hold on
plot3(receiver(:,1), receiver(:,2), receiver(:,3), 'o')
plot3(Q(:,1), Q(:,2), Q(:,3))
plot3(apex(:,1), apex(:,2), apex(:,3), 'o')
title('Scene')
legend('Source', 'Receiver', 'Planes', 'Apex')
view([0 90])
xlim([-4 2])
ylim([-4 2])

nexttile([1 3])
semilogx(fvec, apexTfmag2(:,end))
hold on
semilogx(fvec, planeTfmag2(:,end))
semilogx(fvec, planeApexTfmag2(:,end))
semilogx(fvec, utdTfmag(:,end))
semilogx(fvec, tfmag.diff2)
title('Frequency responses')
xlabel('Frequency')
ylabel('Magnitude')
legend('BTM Apex Daisy Chain', 'BTM Plane Apex Daisy Chain', 'BTM Plane Daisy Chain', 'UTD', 'True BTM', 'Location', 'southwest')
xlim([20 20000])
ylim([-70 0])

%% Test
scale = KimCorrection(data, 2, true);

meanSize = (data.rS + mean(data.W) + data.rR) / 3;
minSize = min(data.rS, min(data.rR, min(data.W)));
maxSize = max(data.rS, max(data.rR, max(data.W)));

