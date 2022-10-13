close all
clear all

fs = 96e3;
nfft = 4096;

barrierRadius = [1, 2];
radiusS = [1, 2];
radiusR = [2, 3];
numEdges = 3;

barrierHeight = 10;
thetaS = 30;
thetaR = 358;
[zS, zR] = deal(barrierHeight / 2);
controlparameters = struct('fs', fs, 'nfft', nfft, 'difforder', 1);

barrierRadius = 1.5;
radiusS = 1.1;
radiusR = 1.4;

[irAprx, tfmagAprx, tvecAprx, fvecAprx, tfcomplexAprx] = SingleNthOrderBarrier1storder(barrierRadius, barrierHeight, thetaS, thetaR, radiusS, radiusR, zS, zR, controlparameters, numEdges, true);

controlparameters = struct('fs', 48e3, 'nfft', nfft, 'difforder', numEdges);

[ir, tfmag, tvec, fvec, tfcomplex] = SingleNthOrderBarrier(barrierRadius, barrierHeight, thetaS, thetaR, radiusS, radiusR, zS, zR, controlparameters, true);

PlotSpectogram([tfcomplexAprx(:,end) tfcomplexAprx(:,end)], fvecAprx, [0 10], [-70 0], 'BTM approach', true, true)

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
semilogx(fvec, tfmag.diff3)
hold on
semilogx(fvecAprx, tfmagAprx(:,end))
for i = 1:numEdges
    semilogx(fvecAprx, tfmagAprx(:,i),'--')
end
hold off
legend('BTM hod', 'BTM 1st order')
xlim([20 20000])
ylim([-100 0])

geometry = GeometryNthOrderBarrier(barrierRadius, radiusS, radiusR);

[result, geometry] = SingleNthOrderBarrierArray(geometry, barrierHeight, thetaS, thetaR, zS, zR, controlparameters);

