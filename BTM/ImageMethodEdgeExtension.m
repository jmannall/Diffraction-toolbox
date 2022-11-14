clear all
close all

wedgeLength = 1;
wedgeIndex = 320;
thetaS = 20;
thetaR = 310;
radiusS = 1;
radiusR = 1;
zS = 0.5;
zR = 0.5;
fs = 96e3;
nfft = 4096;
c = 344;
controlparameters = struct('fs', fs, 'nfft', nfft, 'c', c, 'diffoder', 1, 'saveFiles', 3);

%% Each Path
% Direct
[ir, tfmag, tvec, fvec, tfcomplex] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaR, radiusS, radiusR, zS, zR, controlparameters, false);
tfcomplexAll(:,1) = tfcomplex.diff1;

% sp-ed
zS = -zS;
[ir, tfmag, tvec, fvec, tfcomplex] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaR, radiusS, radiusR, zS, zR, controlparameters, false);
tfcomplexAll(:,2) = tfcomplex.diff1;

% ed-sp
zS = -zS;
zR = -zR;
[ir, tfmag, tvec, fvec, tfcomplex] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaR, radiusS, radiusR, zS, zR, controlparameters, false);
tfcomplexAll(:,3) = tfcomplex.diff1;

% sp-ed-sp
zS = -zS;
[ir, tfmag, tvec, fvec, tfcomplex] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaR, radiusS, radiusR, zS, zR, controlparameters, false);
tfcomplexAll(:,4) = tfcomplex.diff1;

%% Bundled paths
newWedgeLength = 2 * wedgeLength;
zS = -zS + wedgeLength;
zR = -zR + wedgeLength;

% Direct and sp-ed-sp
[ir, tfmag, tvec, fvec, tfcomplex] = SingleWedge(newWedgeLength, wedgeIndex, thetaS, thetaR, radiusS, radiusR, zS, zR, controlparameters, false);
tfcomplexBundle(:,1) = tfcomplex.diff1;

% sp-ed and ed-sp
zS = zS - wedgeLength;
[ir, tfmag, tvec, fvec, tfcomplex] = SingleWedge(newWedgeLength, wedgeIndex, thetaS, thetaR, radiusS, radiusR, zS, zR, controlparameters, false);
tfcomplexBundle(:,2) = tfcomplex.diff1;

%% Figure

disp('Plot')
figure
semilogx(fvec, mag2db(abs(sum(tfcomplexAll, 2))))
hold on
semilogx(fvec, mag2db(abs(sum(tfcomplexBundle, 2))), '--')