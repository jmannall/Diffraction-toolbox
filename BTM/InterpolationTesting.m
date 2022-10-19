close all

fs = 96000;
nfft = 8192;
c = 344;

n = 500;
wedgeLength = 20;
wedgeIndex = 350;
minAngle = 30;
minBendingAngle = 5;

[rS, rR] = deal(0.5);
[zS, zR] = deal(wedgeLength / 2);
controlparameters = struct('nfft', nfft, 'fs', fs, 'difforder', 1);

[result, geometry, pathLength, validPath] = CreateWedgeSweep(wedgeIndex, minAngle, minBendingAngle, rS, rR, controlparameters, n);


bendingAngle = geometry.bendingAngle;
source = geometry.source;
receiver = geometry.receiver;
tfcomplex = [result.tfcomplex];
tfmag = [result.tfmag];
ir = [result.ir];
fvec = result(1).fvec;
limits = [-70 10];

idx = bendingAngle <= 180;
shadowDiff = [tfcomplex.diff1];
shadowDiff(:,idx) = zeros(nfft / 2, sum(idx));

overlap = (max(bendingAngle) - 180) /  4;
validPath.dir = bendingAngle < 180 + overlap;
updateRate = 20;
samplesPerUpdate = fs / updateRate;
audio = ones(samplesPerUpdate * n, 1);
[dir, dirIr] = DelayLine(audio, pathLength.dir, samplesPerUpdate, validPath.dir, c, fs);
%% Audio
close all

% Load audio
[~, ~, filePath] = NormaliseAudio('audio/whiteNoise.wav', nfft);
[audio, fsInput] = LoopAudio(filePath, 13);
audio = 0.8*resample(audio,fs,fsInput);
saveStem = 'whiteNoise';

validPath = bendingAngle > 0;

% Find idx of position before shadow boundary
idx = find(bendingAngle <= 180, 1, 'last');

% Direct and diffracted magnitudes (linear)
tfmagDiff = abs([tfcomplex.diff1]);
tfmagDir = abs([tfcomplex.direct]);

% Direct and diffracted shadow boundary references
SBdiff = tfmagDiff(:, idx + 1);
SBinc = tfmagDir(:,idx);

% Create scaled magnitude response in shadow zone
idx = bendingAngle > 180;
tfmagScaled = zeros(nfft / 2, sum(idx));
for i = 1:nfft / 2
    tfmagScaled(i,:) = SBinc(i) .* tfmagDiff(i,idx) ./ SBdiff(i);
end

% Store correct transfer functions and phase
direct = [tfcomplex.direct];
tfcomplexDiff = [tfcomplex.diff1];
phase = angle(tfcomplexDiff);

% Interpolate between scaled magnitude (linear) and exact response
window = linspace(0, 1, sum(idx));
interpolated = (1 - window) .* tfmagScaled + window .* abs(tfcomplexDiff(:,idx));

% Add phase and combine scaled response and interpolated response with direct component
interpolated = interpolated .* exp(phase(:,idx).*1i);
interpolated = [direct(:,~idx), interpolated];
scaled = tfmagScaled .* exp(phase(:,idx).*1i);
scaled = [direct(:,~idx), scaled];
noInterpolation = [direct(:,~idx), tfcomplexDiff(:,idx)];

PlotSpectrogram(noInterpolation, fvec, bendingAngle, limits, 'SB_No interpolation', false, true, 'Bending Angle')
PlotSpectrogram(interpolated, fvec, bendingAngle, limits, 'SB_Tsingos interpolated response', false, true, 'Bending Angle')
PlotSpectrogram(scaled, fvec, bendingAngle, limits, 'SB_Tsingos scaled response', false, true, 'Bending Angle')
PlotSpectrogram([tfcomplex.complete], fvec, bendingAngle, limits, 'SB_BTM Original response', false, true, 'Bending Angle')
PlotSpectrogram(interpolated ./ [tfcomplex.complete], fvec, bendingAngle, [-10 10], 'SB_BTM vs interpolated error', false, true, 'Bending Angle')

%% Code works up to here to create plots of interpolation
A = abs(interpolated(:,end)) ./ abs(tfcomplexDiff(:,end));

f = rescale(fvec)';
nfilt = 10;
filt = fir2(nfilt,f,A);

irDiff = [ir.diff1];
ndiff = length(irDiff);
irScaled = zeros(ndiff + nfilt, n);
for i = 1:n
    irScaled(:,i) = conv(filt, irDiff(:,i));
end

irDirect = [[ir.direct]; zeros(nfilt, n)];
irContinuous = [irDirect(:,~idx), irScaled(:,idx)];

audioOut = ConvolveIR(audio, irContinuous, validPath, samplesPerUpdate);

audiowrite(['audio\', saveStem, '96k_Continuous.wav'], audioOut, fs)

PlotSpectogramOfWAV(['audio\', saveStem, '96k_Continuous.wav'], limits, nfft)
%% Check point
idx2 = bendingAngle > 180 & bendingAngle < 180 + (max(bendingAngle) - 180) / 2;
window = sqrt([linspace(0, 1, sum(idx2)), ones(1, sum(idx) - sum(idx2))]);
tfmagInterpolate = tfmagDir;
tfmagInterpolate(:,idx) = (1 - window) .* tfmagScaled + window .* tfmagDiff(:,idx);

PlotSpectogram(tfmagInterpolate, fvec, bendingAngle, limits, 'Scaled', false, false, 'Bending Angle')
PlotSpectogram([tfmag.complete], fvec, bendingAngle, limits, 'Original', false, false, 'Bending Angle')
PlotSpectogram([tfmag.complete] - tfmagInterpolate, fvec, bendingAngle, [-5 5], 'Error', false, false, 'Bending Angle')

interpolateTfcomplex = tfmagInterpolate.*exp(angle([tfcomplex.diff1])*1i);
interpolateIr = ifft([interpolateTfcomplex, fliplr(conj(interpolateTfcomplex))]);

diffIr = [ir.diff1];
ndiff = length(diffIr);
ndir = length(dirIr);
idx = bendingAngle <= 180;
diffIr(:,idx) = zeros(ndiff, sum(idx));

padding = zeros(ndiff - ndir, n);
dirIr = [dirIr; padding];

[~, tfcomplexTest] = IrToTf(dirIr + diffIr, nfft);

validPath = bendingAngle > 0;
audioOut = ConvolveIR(audio, [ir.complete], validPath, samplesPerUpdate);
audiowrite(['audio\', saveStem, '96k_Original.wav'], audioOut, fs);

shadowIr = [ir.diff1];
shadowIr(:,idx) = zeros(ndiff, sum(idx));

audioOutNoInt = ConvolveIR(audio, [ir.direct] + shadowIr, validPath, samplesPerUpdate);
audiowrite(['audio\', saveStem, '96k_NoInterpolation.wav'], audioOutNoInt, fs);

idx = bendingAngle > 180 & bendingAngle < 180 + overlap;

interpolateLength = sum(idx);
% window = linspace(0, 1, interpolateLength);
window = hamming(2 * interpolateLength);
window = window(1:interpolateLength)';
dirIr(:,idx) = (1 - window) .* dirIr(:,idx);
diffIr(:,idx) = window .* diffIr(:,idx);

[~, tfcomplexTest2] = IrToTf(dirIr + diffIr, nfft);

audioOutInt = ConvolveIR(audio, dirIr + diffIr, validPath, samplesPerUpdate);

audiowrite(['audio/', saveStem, '96k_Interpolation.wav'], audioOutInt, fs);

%% Plots
close all
nfft = 8192;
PlotSpectogramOfWAV(['audio\', saveStem, '96k_Original.wav'], limits, nfft)
PlotSpectogramOfWAV(['audio\', saveStem, '96k_Interpolation.wav'], limits, nfft)
PlotSpectogramOfWAV(['audio\', saveStem, '96k_NoInterpolation.wav'], limits, nfft)

PlotSpectogram(tfcomplexTest, fvec, bendingAngle, limits, 'Original', false, true, 'Bending Angle')
PlotSpectogram(tfcomplexTest2, fvec, bendingAngle, limits, 'Interpolate', false, true, 'Bending Angle')

PlotSpectogram([tfcomplex.direct] + [tfcomplex.geom] + [tfcomplex.diff1], fvec, bendingAngle, limits, 'Complete', false, true, 'Bending Angle')
% PlotSpectogram([tfcomplex.diff1], fvec, bendingAngle, limits, 'Diffraction', false, false, 'Bending Angle')
% PlotSpectogram(shadowDiff, fvec, bendingAngle, limits, 'Shadow Diffraction', false, false, 'Bending Angle')
% PlotSpectogram([tfcomplex.direct] + [tfcomplex.geom] + shadowDiff, fvec, bendingAngle, limits, 'Complete Shadow Diffraction', false, false, 'Bending Angle')