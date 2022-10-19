close all
clear all

set(0, 'DefaultLineLineWidth', 2);

fs = 96e3;
nfft = 4096;

barrierRadius = [1, 2];
radiusS = [1, 2];
radiusR = [2, 3];
numEdges = 2;

barrierHeight = 10;
thetaS = 30;
thetaR = 358;
[zS, zR] = deal(barrierHeight / 2);
controlparameters = struct('fs', fs, 'nfft', nfft, 'difforder', 1);

barrierRadius = 1.5;
radiusS = 1.5;
radiusR = 1.1;

thetaS = 45;
thetaR = 205;
barrierRadius = 0.89 / 2;
radiusS = 0.8;
radiusR = 0.23;

[irAprx, tfmagAprx, tvecAprx, fvecAprx, tfcomplexAprx] = SingleNthOrderBarrier1storder(barrierRadius, barrierHeight, thetaS, thetaR, radiusS, radiusR, zS, zR, controlparameters, numEdges, true);

controlparameters = struct('fs', fs, 'nfft', nfft, 'difforder', numEdges);

[ir, tfmag, tvec, fvec, tfcomplex] = SingleNthOrderBarrier(barrierRadius, barrierHeight, thetaS, thetaR, radiusS, radiusR, zS, zR, controlparameters, true);

PlotSpectrogram([tfcomplexAprx(:,end) tfcomplexAprx(:,end)], fvecAprx, [0 10], [-70 0], 'BTM approach', true, true)

% figure
% semilogx(fvec, tfmag.diff2)
% hold on
% semilogx(fvec1, tfmagAprx)
% semilogx(fvec1, tfmag1.diff1,'--')
% semilogx(fvec2, tfmag2.diff1,'--')
% hold off
% legend('BTM hod', 'BTM 1st order')
% xlim([20 20000])
% ylim([-100 0])    

figure
semilogx(fvec, tfmag.diff2)
hold on
semilogx(fvecAprx, tfmagAprx(:,end))
for i = 1:numEdges
    semilogx(fvecAprx, tfmagAprx(:,i),'--')
end
hold off
legend('BTM hod', 'BTM 1st order')
xlim([20 20000])
ylim([-100 0])



