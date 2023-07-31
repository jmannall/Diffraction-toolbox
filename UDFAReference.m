
close all

fs = 48e3;
nfft = 8192;
c = 344;

controlparameters = struct('fs', fs, 'nfft', nfft', 'c', c, 'noDirect', false, 'saveFiles', 0);

wedgeLength = 4;
wedgeIndex = 270;
thetaS = 10;
thetaR = 250;
rS = 2;
rR = 1;
zS = 1;
zR = 1;
n = 4;

[tfmag, fvec] = SingleIIRWedge(wedgeLength, wedgeIndex, thetaS, thetaR, rS, rR, zS, zR, n, controlparameters);

controlparameters.fs = 2 * fs;
controlparameters.nfft = 2 * nfft;
[~, tfmagBtm, ~, fvecBtm] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaR, rS, rR, zS, zR, controlparameters, false);
%%
input = extractdata(NNInputFromGeometry(wedgeIndex, wedgeLength, thetaR, thetaS, rS, rR, zS, zR, false))';
dZ = abs(zR - zS);
pathLength = sqrt((rS + rR) .^ 2 + dZ .^ 2);

[bestz, bestp, bestk] = myBestNN(input);
[tfmagNNBest, ~] = CreateIIRFilter(bestz, bestp, bestk, nfft, fs);
tfmagNNBest = mag2db(db2mag(tfmagNNBest) / pathLength);


[smallz, smallp, smallk] = mySmallNN(input);
[tfmagNNSmall, ~] = CreateIIRFilter(smallz, smallp, smallk, nfft, fs);
[tfmagNNSmall_Test, ~] = CreateIIRFilter([0.905099750; 0.146374494], [0.957741261; 0.640204012], 0.100753546, nfft, fs);
tfmagNNSmall_Test = mag2db(db2mag(tfmagNNSmall_Test) / pathLength);
[b, a] = zp2tf(smallz, smallp, smallk);
tfmagNNSmall = mag2db(db2mag(tfmagNNSmall) / pathLength);

%%
fc = [176, 775, 3408];
for i = 1:4
    if i == 1
        fMin = 20;
    else
        fMin = fc(i - 1);
    end
    if i == 4
        fMax = 20e3;
    else
        fMax = fc(i);
    end
    midPoints(i) = sqrt(fMin .* fMax);
end

[~, phii] = CalculateApex(rS, rR, zS, zR, wedgeLength, true);
controlparameters.fvec = midPoints;
[tfmagUtd, fvecUtd, tfcomplexUtd] = SingleUTDWedge(thetaS, thetaR, rS, rR, wedgeIndex, phii, controlparameters); % 174
[tfmagUtd, fvecUtd, tfcomplexUtd] = SingleUTDWedgeInterpolated(thetaS, thetaR, rS, rR, wedgeIndex, phii, controlparameters); % 174

tfUtd = db2mag(tfmagUtd);

LinkwitzRiley(fvec, fs, nfft, tfmagUtd);
%%
figure
semilogx(fvecBtm, tfmagBtm.diff1)
hold on
%semilogx(fvec, tfmag, '--')
%semilogx(fvec, tfmagNNBest, ':')
semilogx(fvec, tfmagNNSmall, '--')
semilogx(fvec, tfmagNNSmall_Test, ':')
semilogx(fvecUtd, tfmagUtd, '--')
grid on
xlim([20 20e3])