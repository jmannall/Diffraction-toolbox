close all
clear all

fileStem = 'audio\whiteNoise';
audioFilePath  = [fileStem, '.wav'];
exampleLength = 10;
nfft = 8192;
fs = 48e3;
c = 344;    % speed of sound

% [audio, fs] = NormaliseAudio(audioFilePath, nfft);
% 
% audiowrite([fileStem, '_Normalised.wav'], audio, fs)
% 
% [audio, fs, audioLength] = LoopAudio([fileStem, '_Normalised.wav'], exampleLength);
% 
% updateRate = 20;
% samplesPerUpdate = fs / updateRate;
% 
% numWindows = 2 * floor(audioLength / samplesPerUpdate) - 2;

wedgeIndex = 350;
minBendingAngle = 160;
minAngle = 20;
rS = 1;
rR = 1;
numWindows = wedgeIndex - (minAngle + minBendingAngle);

%% BTM

controlparameters = struct('nfft', 2 * nfft, 'fs', 2 * fs, 'c', c, 'diffOrder', 1, 'saveFiles', 2, 'noDirect', false);
[result, geometry, pathLength, validPath] = CreateWedgeSweep(wedgeIndex, minAngle, minBendingAngle, rS, rR, controlparameters, numWindows, false);

ir = [result.ir];
tfmag = [result.tfmag];
tfcomplex = [result.tfcomplex];
fvec = result(1).fvec;
bendingAngle = geometry.bendingAngle';

fvec2 = fs/nfft*[0:nfft/2-1];

%% BTM interpolated
[result, geometry, ~, ~] = CreateWedgeSweep(wedgeIndex, minAngle, minBendingAngle, rS, rR, controlparameters, numWindows, true);
irI = [result.ir];
tfmagI = [result.tfmag];
tfcomplexI = [result.tfcomplex];

tfcomplexI(:,validPath.dir) = 1 ./ pathLength.dir(validPath.dir) .* tfcomplexI(:,validPath.dir);
tfcomplexI(:,validPath.diffShadow) = 1 ./ pathLength.diff(validPath.diffShadow) .* tfcomplexI(:,validPath.diffShadow);

%% Neural network

loadPath = ['NNSaves', filesep, 'iir-2057_0001-1-09-099-3-25.mat'];

load(loadPath, "net");

wedgeLength = 20;
[zS, zR] = deal(wedgeLength / 2);
rS = 1;
rR = 1;
thetaS = minAngle;

const = ones(size(geometry.receiver));
rR = const .* rR;
zR = const .* zR;

output.NN = CreateNNOutput(net, wedgeIndex, wedgeLength, geometry.receiver, thetaS, rS, rR, zS, zR);

biquad = false;
[tfcomplexNN, b.NN, a.NN] = CalculateNN(output.NN, [tfcomplex.direct], validPath.diffShadow, pathLength.diff', nfft, fs, biquad);

%[audioPath.NN, ir.NN] = DelayLineIIRFilter(audio, [pathLength.diff, pathLength.dir], windowLength, validPath.diff, b.NN, a.NN, c, fs, validPath.diffShadow);

%[~, tfcomplex.NN] = IrToTf(ir.NN, nfft);


%numBiquads = 2;
%NNidx = '015b';
%[tfmagNN, tfcomplexNN, b, a,] = ProcessBiquadNNOutputWedgeSweep(NNidx, wedgeIndex, bendingAngle, minAngle, numBiquads, pathLength.diff, c, nfft, fs);

%dir = DelayLine(audio, pathLength.dir, samplesPerUpdate, validPath.dir, c, fs);
%spec = DelayLine(audio, pathLength.spec, samplesPerUpdate, validPath.spec, c, fs);
%diffNN = DelayLineBiquad(audio, pathLength.diff, samplesPerUpdate, validPath.diff, b, a, numBiquads, c, fs);
%BTM = ConvolveIR(audio, [ir.diff1], samplesPerUpdate, validPath.diff);
% BTMcomplete= ConvolveIR(audio, [ir.complete], samplesPerUpdate, validPath.diff);

%disp('Write audio files')
%audiowrite([fileStem, '_dir.wav'], dir, fs)
%audiowrite([fileStem, '_spec.wav'], spec, fs)
%audiowrite([fileStem, '_NN', NNidx, 'diff.wav'], diffNN, fs)
%audiowrite([fileStem, '_NN', NNidx, 'All.wav'], dir + spec + diffNN, fs)
%audiowrite([fileStem, '_BTM.wav'], BTM, fs)
%audiowrite([fileStem, '_BTMAll.wav'], dir + spec + BTM, fs)
% audiowrite([fileStem, '_BTMTest.wav'], BTMcomplete, fs)

%% Plots
close all

disp('Create plots')
limits = [-70 0];
% PlotSpectogramOfWAV([fileStem, '_dir.wav'], limits, nfft)
% PlotSpectogramOfWAV([fileStem, '_spec.wav'], limits, nfft)
% PlotSpectogramOfWAV([fileStem, '_NNdiff.wav'], limits, nfft)
% PlotSpectogramOfWAV([fileStem, '_NNAll.wav'], limits, nfft)
% PlotSpectogramOfWAV([fileStem, '_BTM.wav'], limits, nfft)
% PlotSpectogramOfWAV([fileStem, '_BTMAll.wav'], limits, nfft)
% PlotSpectogramOfWAV([fileStem, '_BTMTest.wav'], limits, nfft)

phase = false;
save = true;
PlotSpectrogram([tfcomplex.direct], fvec, bendingAngle, limits, 'Direct', phase, save, 'Bending Angle')
PlotSpectrogram([tfcomplex.geom], fvec, bendingAngle, limits, 'Specular', phase, save, 'Bending Angle')
%PlotSpectogram([tfcomplex.direct] + [tfcomplex.geom] + tfcomplexNN, fvec, bendingAngle, limits, ['NN', NNidx, ' All'], false, 'Bending Angle')
PlotSpectrogram([tfcomplex.diff1], fvec, bendingAngle, limits, 'BTM diffraction', phase, save, 'Bending Angle')
PlotSpectrogram([tfcomplex.complete], fvec, bendingAngle, limits, 'BTM All', phase, save, 'Bending Angle')
PlotSpectrogram(tfcomplexI + [tfcomplex.geom], fvec, bendingAngle, limits, 'BTM Interpolated', phase, save, 'Bending Angle')
geomNN = [tfcomplex.geom];
PlotSpectrogram(tfcomplexNN + geomNN(1:end/2,:), fvec(1:end/2), bendingAngle, limits, 'NN diffraction', phase, save, 'Bending Angle')

PlotSpectrogram(tfcomplexNN ./ tfcomplexI(1:end/2,:), fvec(1:end/2), bendingAngle, [-5 5], 'NN error', phase, save, 'Bending Angle')