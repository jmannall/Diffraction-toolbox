clear all
close all

timeStep = 0.08;   % Length of time spent at each position

zR = [0.2, 0.1; 0.5, 0.4];
zI = [0.01, 0.01; 0.01, 0.01];
pR = [0.8, 0.7; 0.1, 0];
pI = [0.01, 0.01; 0.02, 0.02];
k = [0.2, 0.1];

pathLengthDiff = [2, 2];
%pathLengthDir = 2:-0.01:1;

numFreq = 4096;

numBiquads = size(zR,1);

[audio, fs] = audioread("sineTone.wav");
step = round(timeStep * fs);   % Number of samples at each position
updateRate = 20;
bufferLength = fs / updateRate;
n = updateRate / 2;
pathLengthDir = linspace(1, 3, n);
numObservations = size(zR,2);
numObservations = size(pathLengthDir,2);
impulseLength = numObservations * bufferLength;

c = 343;    % speed of sound
[delayLengthDiff, fracDelayDiff] = CalculateDelay(pathLengthDiff, c, fs);
[delayLengthDir, fracDelayDir] = CalculateDelay(pathLengthDir, c, fs);

% [tfmag, fvec] = CreateBiquad(zR, zI, pR, pI, k, numFreq, fs);
% 
% figure
% semilogx(fvec, extractdata(tfmag))

% [ir, tfmag, tfcomplex] = ComputeWedgeResponses(audio, 8192);

% figure
% semilogx(fvec, tfmag)

figure
plot(1:impulseLength, audio(1:impulseLength))

% [b, a] = deal(zeros(3,numBiquads,numObservations));
% for i = 1:numBiquads
%     b(:,i,:) = [ones(1,numObservations); -2 * zR(i,:); zR(i,:).^2 + zI(i,:).^2];
%     a(:,i,:) = [ones(1,numObservations); -2 * pR(i,:); pR(i,:).^2 + pI(i,:).^2];
% end
% b(:,1,:) = squeeze(b(:,1,:)) .* k;

inputBuffer = zeros(3, numBiquads);     % Diffraction IIR filter input buffer
outputBuffer = zeros(3, numBiquads);      % Diffraction IIR filter output buffer



fracFilterLength = 32;

[delayLengthDir, fracDelayDir] = CalculateDelay(pathLengthDir, c, fs);

bufferLengthDir = max(delayLengthDir);
bufferDir = zeros(bufferLengthDir, 1);
readDir = delayLengthDir;
pathLength = pathLengthDir(1);
writeFilter = 0;
readFilter = 31:-1:0;
bufferFilter = zeros(fracFilterLength, 1);
overlap = floor(bufferLength / 2);
window = hamming(2 * bufferLength);

x = audio;      % input
y = zeros(1, numObservations * bufferLength);   % output

writeDir = 0;
[delayLengthDir1, fracDelayDir1] = CalculateDelay(pathLengthDir(2), c, fs);
[readDir1, readDir2] = deal(mod(delayLengthDir1, bufferLengthDir));

for k = 2:numObservations
    delayLengthDir2 = delayLengthDir1;
    fracDelayDir2 = fracDelayDir1;
    [delayLengthDir1, fracDelayDir1] = CalculateDelay(pathLengthDir(k), c, fs);
    v = delayLengthDir2 - delayLengthDir1;
    readDir2 = readDir1 + v;
    for i = 1:bufferLength
        idx = (k - 1) * bufferLength + i;
        
        readDir2 = mod(readDir2, bufferLengthDir) + 1;
        readDir1 = mod(readDir1, bufferLengthDir) + 1;
        writeDir = mod(writeDir, bufferLengthDir) + 1;
        y(idx) = window(bufferLength + i) * bufferDir(readDir1) + window(i) * bufferDir(readDir2);
        bufferDir(writeDir) = x(idx);
    end
end

audiowrite("sineTone_filt.wav", y(1:impulseLength), fs)

figure
plot(1:impulseLength, y(1:impulseLength))




i0var = 1;
FDs = i0var+5*(0:0.2:0.8);
for k = 2:numObservations
%     diffBuffer = zeros(sampleDelayDiff(k), 1);
    delayLengthDirOld = delayLengthDir;
    [delayLengthDir, fracDelayDir] = CalculateDelay(pathLengthDir(k), c, fs);
    [hDir, iDir, bwDir] = designFracDelayFIR(fracDelayDir, fracFilterLength);
    v = delayLengthDirOld - delayLengthDir;
    readDir = readDir + v - overlap;
%     v = (pathLengthDir(k) - pathLengthDir(max(k - 1, 1))) / timeStep;
%     g = -v / c;
    idx = max((k - 1) * (bufferLength - overlap) + 1, 1):k * (bufferLength - overlap) + overlap;
    input = x(idx); 
    for i = 1:bufferLength
        idx = (k - 1) * (bufferLength - overlap) + i;   % Input sample idx
        writeDir = mod(idx - 1, bufferLengthDir) + 1;
        bufferDir(writeDir) = bufferDir(writeDir) + window(i) * x(idx);
x
        readDir = mod(readDir, bufferLengthDir) + 1;
        writeDir = mod(idx - 1, bufferLengthDir) + 1;

        readFilter = mod(readFilter, fracFilterLength) + 1;
        writeFilter = mod(writeFilter, fracFilterLength) + 1;

        bufferFilter(writeFilter) = bufferDir(readDir);

        dir = (1 / pathLengthDir(k)) * hDir * bufferFilter(readFilter);
        bufferDir(writeDir) = x(idx);

        %dirBuffer = circshift(dirBuffer, 1);
%         diffBuffer = circshift(diffBuffer, 1);
%         inputBuffer = circshift(inputBuffer, 1);
%         outputBuffer = circshift(outputBuffer, 1);
%         pathLength = pathLength - g;
%         [delayDir, fracDelayDir] = CalculateDelay(pathLength, c, fs);
        %readDir = mod(readDir, bufferLengthDir) + 1 + g;
%         [hDir, iDir, bwDir] = designFracDelayFIR(fracDelayDir,fracFilterLength);
%         readDir = (delayDir:delayDir + fracFilterLength - 1) - iDir + 1;
%         inputBuffer = bufferDir(mod(readDir, bufferLengthDir) + 1);
%         for j = 1:fracFilterLength
%             output = inputBuffer .* hDir;
%         end

        %scale = readDir - readDirMin;
        %readDirMin = mod(readDirMin - 1, bufferLengthDir) + 1;
        %dir = (1 / pathLengthDir(k)) * ((1 - scale) * bufferDir(readDirMin) + scale * bufferDir(readDirMax));
        %dir = (1 / pathLengthDir(k)) * bufferDir(mod(round(readDir), bufferLengthDir) + 1);
%         diffBuffer(1,1) = x(idx);
    
%         diff = diffBuffer(sampleDelayDiff(k)); 
        % Update the input and output buffer for each biquad
%         for j = 1:numBiquads
%             inputBuffer(1,j) = diff;
%             diff = b(1,j,k) * inputBuffer(1,j) + b(2,j,k) * inputBuffer(2,j) + b(3,j,k) * inputBuffer(3,j) - a(2,j,k) * outputBuffer(2,j) - a(3,j,k) * outputBuffer(3,j);
%             outputBuffer(1,j) = diff;
%         end
        %bufferDir(writeDir) = x(idx);
        y(idx) = y(idx) + window(i) * dir;
    end
end

[ir, tfmag, tfcomplex] = ComputeWedgeResponses(y', 8192);

% figure
% semilogx(fvec, tfmag)
 
figure
plot(1:impulseLength, y(1:impulseLength))

audiowrite("sineTone_filt.wav", y(1:impulseLength), fs)
