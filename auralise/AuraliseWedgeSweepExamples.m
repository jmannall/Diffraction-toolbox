close all
clear all

fileStem = 'audio\whiteNoise';
audioFilePath  = [fileStem, '.wav'];
exampleLength = 10;
nfft = 8192;

[audio, fs] = NormaliseAudio(audioFilePath, nfft);

audiowrite([fileStem, '_Normalised.wav'], audio, fs)

[audio, fs, audioLength] = LoopAudio([fileStem, '_Normalised.wav'], exampleLength);

updateRate = 20;
samplesPerUpdate = fs / updateRate;

numWindows = 2 * floor(audioLength / samplesPerUpdate) - 2;

wedgeIndex = 350;
minBendingAngle = 5;
minAngle = 30;
rS = 1;
rR = 1;

controlparameters = struct('nfft', nfft, 'fs', fs);
[result, geometry, pathLength, validPath] = CreateWedgeSweep(wedgeIndex, minAngle, minBendingAngle, rS, rR, controlparameters, numWindows);

ir = [result.ir];
tfmag = [result.tfmag];
tfcomplex = [result.tfcomplex];
fvec = result(1).fvec;
bendingAngle = geometry.bendingAngle';

c = 344;    % speed of sound

numBiquads = 2;
NNidx = '015b';
[tfmagNN, tfcomplexNN, b, a,] = ProcessBiquadNNOutputWedgeSweep(NNidx, wedgeIndex, bendingAngle, minAngle, numBiquads, pathLength.diff, c, nfft, fs);

dir = DelayLine(audio, pathLength.dir, samplesPerUpdate, validPath.dir, c, fs);
spec = DelayLine(audio, pathLength.spec, samplesPerUpdate, validPath.spec, c, fs);
diffNN = DelayLineBiquad(audio, pathLength.diff, samplesPerUpdate, validPath.diff, b, a, numBiquads, c, fs);
BTM = ConvolveIR(audio, [ir.diff1], samplesPerUpdate, validPath.diff);
% BTMcomplete= ConvolveIR(audio, [ir.complete], samplesPerUpdate, validPath.diff);

disp('Write audio files')
audiowrite([fileStem, '_dir.wav'], dir, fs)
audiowrite([fileStem, '_spec.wav'], spec, fs)
audiowrite([fileStem, '_NN', NNidx, 'diff.wav'], diffNN, fs)
audiowrite([fileStem, '_NN', NNidx, 'All.wav'], dir + spec + diffNN, fs)
audiowrite([fileStem, '_BTM.wav'], BTM, fs)
audiowrite([fileStem, '_BTMAll.wav'], dir + spec + BTM, fs)
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

PlotSpectogram([tfcomplex.direct], fvec, bendingAngle, limits, 'Direct', false, 'Bending Angle')
PlotSpectogram([tfcomplex.geom], fvec, bendingAngle, limits, 'Specular', false, 'Bending Angle')
PlotSpectogram(tfcomplexNN, fvec, bendingAngle, limits, ['NN', NNidx, ' diffraction'], false, 'Bending Angle')
PlotSpectogram([tfcomplex.direct] + [tfcomplex.geom] + tfcomplexNN, fvec, bendingAngle, limits, ['NN', NNidx, ' All'], false, 'Bending Angle')
PlotSpectogram([tfcomplex.diff1], fvec, bendingAngle, limits, 'BTM diffraction', false, 'Bending Angle')
PlotSpectogram([tfcomplex.complete], fvec, bendingAngle, limits, 'BTM All', false, 'Bending Angle')
