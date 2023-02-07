wedgeIndex = 320;

thetaR = 300;
thetaS = 1e-3;

wedgeLength = 20;
rS = 1;
rR = 1;
zS = 5;
zR = 5;

[tfmag, fvec] = SingleWedgeInterpolated(wedgeLength, wedgeIndex, thetaS, thetaR, rS, rR, zS, zR, controlparameters, false);

thetaR = 320 - 1e-3;
thetaS = 20;

tfmag2 = SingleWedgeInterpolated(wedgeLength, wedgeIndex, thetaS, thetaR, rS, rR, zS, zR, controlparameters, false);

figure
semilogx(fvec, tfmag)
hold on
semilogx(fvec, tfmag2, '--')
legend('NN learnt', 'What we want')