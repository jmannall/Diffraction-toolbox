close all
clear all

audioFilePath  = 'whiteNoise';
[audio, fs] = audioread([audioFilePath, '.wav']);
[audio2, ~] = audioread([audioFilePath, '.wav']);
audio = [audio; audio2];

audio = audio * 0.15;
updateRate = 50;
samplesPerUpdate = fs / updateRate;
overlap = floor(samplesPerUpdate / 2);
window = hanning(samplesPerUpdate);

n = floor(2 * length(audio) / samplesPerUpdate - 2);
%n = 798;
rS = 1;
rR = 1;

wedgeIndex = 350;
minAngle = 30.1;
bendingAngle = linspace(4.99, wedgeIndex - minAngle - 0.01, n);

geometry = GeometryAll(wedgeIndex, bendingAngle, minAngle);

disp('Run')
[result, geometry] = SingleWedgeArray(geometry, 20,1,1,10,10,fs);

tfmag = [result.tfmag];
tfcomplex = [result.tfcomplex];
fvec = result(1).fvec;
bendingAngle = geometry.bendingAngle';

c = 344;    % speed of sound

source = [rS * sind(minAngle); rS * cosd(minAngle)];
receiver = [rR * sind(minAngle + bendingAngle); rR * cosd(minAngle + bendingAngle)];

pathLengthDir = vecnorm(source - receiver);
idx = bendingAngle > 180;
pathLengthDir(idx) = 0.001;
pathLengthSpec = vecnorm([-source(1); source(2)] - receiver);
idx = bendingAngle > 180 - 2 * minAngle;
pathLengthSpec(idx) = 0.001;
pathLengthDiff = linspace(rS + rR, rS + rR, n);

[delayDir, fracDelayDir, bufferLengthDir, readDir, writeDir, bufferDir] = CreateSoundComponent(pathLengthDir, c, fs, samplesPerUpdate);
[delaySpec, fracDelaySpec, bufferLengthSpec, readSpec, writeSpec, bufferSpec] = CreateSoundComponent(pathLengthSpec, c, fs, samplesPerUpdate);
[delayDiff, fracDelayDiff, bufferLengthDiff, readDiff, writeDiff, bufferDiff] = CreateSoundComponent(pathLengthDiff, c, fs, samplesPerUpdate);

% Create output variables
[outNN, outNoWindow, outBTM, outNoDiff, outBTMComplete] = deal(zeros(length(audio), 1));
numBuffers = n;

load("NeuralNetwork_015b.mat", 'net')

numBiquads = 2;
numFreq = 4096;

const = ones(1, n);
output = predict(net, dlarray([const * wedgeIndex; bendingAngle; const * minAngle], "CB"));

[zR, zI, pR, pI, k] = CreateZPKFromNNOutput(output, numBiquads);

idx = bendingAngle < 180 - 2 * minAngle | (180 - minAngle < bendingAngle & bendingAngle < 180);
k(idx) = -k(idx);

% [tfmagNN, fvecNN, tfcomplexNN] = CreateBiquad(zR, zI, pR, pI, k, numFreq, fs);

% z = double([extractdata(zR) + 1i*extractdata(zI)]);
% p = double([extractdata(pR) + 1i*extractdata(pI)]);
% k = double(extractdata(k));

biquadSize = 3;
[b, a] = deal(zeros(biquadSize,numBiquads,n));
for i = 1:numBiquads
    b(:,i,:) = [ones(1,n); -2 * zR(i,:); zR(i,:).^2 + zI(i,:).^2];
    a(:,i,:) = [ones(1,n); -2 * pR(i,:); pR(i,:).^2 + pI(i,:).^2];
end
b(:,1,:) = squeeze(b(:,1,:)) .* k;

inputBufferIIR = zeros(biquadSize, numBiquads);     % Diffraction IIR filter input buffer
outputBufferIIR = zeros(biquadSize, numBiquads);      % Diffraction IIR filter output buffer

[inputBufferDir, inputBufferSpec] = deal(zeros(2, 1));

ir = [result.ir];
irBTM = [ir.diff];
irDirSpec = [ir.direct] + [ir.geom];
irBTMComplete = [ir.complete];

bufferLengthBTM = size(irBTM, 1);
bufferBTM = zeros(bufferLengthBTM, 1);
writeBTM = 0;
readBTM = bufferLengthBTM - 1:-1:0;

disp('Loop through buffers')

% fracFilterLength = 32;
% 
% delay = pathLengthDir / c;     % in seconds // time = distance / speed
% sampleDelay = delay * fs;
% fracDelayDir = rem(sampleDelay, 1);
% 
% delay = pathLengthSpec / c;     % in seconds // time = distance / speed
% sampleDelay = delay * fs;
% fracDelaySpec = rem(sampleDelay, 1);
% 
% delay = pathLengthDiff / c;     % in seconds // time = distance / speed
% sampleDelay = delay * fs;
% fracDelayDiff = rem(sampleDelay, 1);

% [fracInputBufferDir, fracInputBufferSpec, fracInputBufferDiff] = deal(zeros(1,fracFilterLength));     % Fractional delay FIR filter input buffer

% readBTM = readBTM - fracFilterLength / 2 + 1;

irDirTest = zeros(bufferLengthBTM, numBuffers);
for i = 1:n
    irDirTest(delayDir(i):delayDir(i)+1,i) = (1 / pathLengthDir(i)) * [(1 - fracDelayDir(i)), fracDelayDir(i)];
end

for k = 1:numBuffers
    [readDir, writeDir] = UpdateReadWrite(readDir, writeDir, overlap, delayDir, k);
    [readSpec, writeSpec] = UpdateReadWrite(readSpec, writeSpec, overlap, delaySpec, k);
    [readDiff, writeDiff] = UpdateReadWrite(readDiff, writeDiff, overlap, delayDiff, k);
%     inputBuffer = zeros(biquadSize, numBiquads);     % Diffraction IIR filter input buffer
%     outputBuffer = zeros(biquadSize, numBiquads);      % Diffraction IIR filter output buffer
    readBTM = readBTM - overlap;
    writeBTM = writeBTM - overlap;
    for i = 1:samplesPerUpdate
%         hDir = designFracDelayFIR(fracDelayDir(k), fracFilterLength);
%         hSpec = designFracDelayFIR(fracDelaySpec(k), fracFilterLength);
%         hDiff = designFracDelayFIR(fracDelayDiff(k), fracFilterLength);
        idx = (k - 1) * (samplesPerUpdate - overlap) + i;
        [readDir, writeDir] = IncrementReadWrite(readDir, writeDir, bufferLengthDir);
        [readSpec, writeSpec] = IncrementReadWrite(readSpec, writeSpec, bufferLengthSpec);
        [readDiff, writeDiff] = IncrementReadWrite(readDiff, writeDiff, bufferLengthDiff);
        [readBTM, writeBTM] = IncrementReadWrite(readBTM, writeBTM, bufferLengthBTM);
        inputBufferIIR = circshift(inputBufferIIR, 1);
        outputBufferIIR = circshift(outputBufferIIR, 1);
        inputBufferDir = circshift(inputBufferDir, 1);
        inputBufferSpec = circshift(inputBufferSpec, 1);
%         fracInputBufferDir = circshift(fracInputBufferDir, 1);
%         fracInputBufferSpec = circshift(fracInputBufferSpec, 1);
%         fracInputBufferDiff = circshift(fracInputBufferDiff, 1);
        if bendingAngle(k) <= 180
            dir = (1 / pathLengthDir(k)) * bufferDir(readDir);
        else
            dir = 0;
        end
        if bendingAngle(k) <= (180 - 2 * minAngle)
            spec = (1 / pathLengthSpec(k)) * bufferSpec(readSpec);
        else
            spec = 0;
        end
        inputBufferDir(1) = dir;
        inputBufferSpec(1) = spec;
        dir = inputBufferDir(1) * (1 - fracDelayDir(k)) + inputBufferDir(2) * fracDelayDir(k); % A(1-r)*d(n-n0)+A*r*d(n-n0-1);
        spec = inputBufferSpec(1) * (1 - fracDelaySpec(k)) + inputBufferSpec(2) * fracDelaySpec(k); % A(1-r)*d(n-n0)+A*r*d(n-n0-1);

        diff = (1 / pathLengthDiff(k)) * bufferDiff(readDiff);
%         for j = 1:fracFilterLength
%             fracInputBufferDir(1,j) = dir;
%             dirFrac = sum(hDir .* fracInputBufferDir);
%         end
%         for j = 1:fracFilterLength
%             fracInputBufferSpec(1,j) = spec;
%             specFrac = sum(hSpec .* fracInputBufferSpec);
%         end
%         for j = 1:fracFilterLength
%             fracInputBufferDiff(1,j) = diff;
%             diffFrac = sum(hDiff .* fracInputBufferDiff);
%         end
        for j = 1:numBiquads
            inputBufferIIR(1,j) = diff;
            diff = b(1,j,k) * inputBufferIIR(1,j) + b(2,j,k) * inputBufferIIR(2,j) + b(3,j,k) * inputBufferIIR(3,j) - a(2,j,k) * outputBufferIIR(2,j) - a(3,j,k) * outputBufferIIR(3,j);
            outputBufferIIR(1,j) = diff;
        end
        BTM = bufferBTM(readBTM)' * irBTM(:, k);
        if bendingAngle(k) <= 180
            diffShadow = 0;
        else
            diffShadow = BTM;
        end
        dirSpec = bufferBTM(readBTM)' * irDirSpec(:,k);
        BTMComplete = bufferBTM(readBTM)' * irBTMComplete(:, k);
        outNN(idx) = outNN(idx) + window(i) * (dir + spec + diff);
        outBTM(idx) = outBTM(idx) + window(i) * (dir + spec + BTM);
        outNoDiff(idx) = outNoDiff(idx) + window(i) * (dirSpec + diffShadow);
        outBTMComplete(idx) = outBTMComplete(idx) + window(i) * (dirSpec + BTM);
        bufferDir(writeDir) = audio(idx);
        bufferSpec(writeSpec) = audio(idx);
        bufferDiff(writeDiff) = audio(idx);
        bufferBTM(writeBTM) = audio(idx);
    end
end

%% Plots

close all
disp('Plot')

x = 1:length(outNN);
compare = audio(x);
figure
plot(x, [compare outBTMComplete outBTM outNN])
legend('Original', 'BTM Complete', 'BTM', 'NN', 'Location', 'northoutside')
xlim([0 500])

audiowrite([audioFilePath, '_NN.wav'], outNN, fs);
audiowrite([audioFilePath, '_BTM.wav'], outBTM, fs);
audiowrite([audioFilePath, '_NoDiff.wav'], outNoDiff, fs);
audiowrite([audioFilePath, '_BTMComplete.wav'], outBTMComplete, fs);

sg = 400;
ov = 300;

limits = [-70 0];

[x, fs] = audioread([audioFilePath, '.wav']);
[S,F,T] = spectrogram(x,sg,ov,[],fs,"yaxis");
figure
sh = surf(T,F,20*log10(abs(S)));
colorbar
view([0 90])
axis tight
xlabel('Time')
ylabel('Frequency')
set(gca,'YScale','log')
set(sh, 'EdgeColor','none')
clim(limits)
title('Original')

[x, fs] = audioread([audioFilePath, '_NN.wav']);
[S,F,T] = spectrogram(x,sg,ov,[],fs,"yaxis");
figure
sh = surf(T,F,20*log10(abs(S)));
colorbar
view([0 90])
axis tight
xlabel('Time')
ylabel('Frequency')
set(gca,'YScale','log')
set(sh, 'EdgeColor','none')
clim(limits)
title('Neural Network')

[x, fs] = audioread([audioFilePath, '_BTM.wav']);
[S,F,T] = spectrogram(x,sg,ov,[],fs,"yaxis");
figure
sh = surf(T,F,20*log10(abs(S)));
colorbar
view([0 90])
axis tight
grid on
xlabel('Time')
ylabel('Frequency')
set(gca,'YScale','log')
set(sh, 'EdgeColor','none')
clim(limits)
title('BTM')

[x, fs] = audioread([audioFilePath, '_NoDiff.wav']);
[S,F,T] = spectrogram(x,sg,ov,[],fs,"yaxis");
figure
sh = surf(T,F,20*log10(abs(S)));
colorbar
view([0 90])
axis tight
grid on
xlabel('Time')
ylabel('Frequency')
set(gca,'YScale','log')
set(sh, 'EdgeColor','none')
clim(limits)
title('Shadow only Diffraction')

[x, fs] = audioread([audioFilePath, '_BTMComplete.wav']);
[S,F,T] = spectrogram(x,sg,ov,[],fs,"yaxis");
figure
sh = surf(T,F,20*log10(abs(S)));
colorbar
view([0 90])
axis tight
grid on
xlabel('Time')
ylabel('Frequency')
set(gca,'YScale','log')
set(sh, 'EdgeColor','none')
clim(limits)
title('BTM Complete')
