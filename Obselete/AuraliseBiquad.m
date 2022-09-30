clear all
close all

timeStep = 1;   % Length of time spent at each position

zR = [0.2, 0.1; 0.5, 0.4];
zI = [0.01, 0.01; 0.01, 0.01];
pR = [0.8, 0.7; 0.1, 0];
pI = [0.01, 0.01; 0.02, 0.02];
k = [0.2, 0.1];

pathLength = 2;

step = round(timeStep * fs);   % Number of samples at each position
numFreq = 4096;
numObservations = size(zR,2);
numBiquads = size(zR,1);

[audio, fs] = audioread("whiteNoise.wav");
impulseLength = min(numObservations,length(audio));

c = 343;    % speed of sound
delay = pathLength / c;     % in seconds // time = distance / speed
sampleDelay = delay * fs;

delayLength

% [tfmag, fvec] = CreateBiquad(zR, zI, pR, pI, k, numFreq, fs);
% 
% figure
% semilogx(fvec, extractdata(tfmag))

% [ir, tfmag, tfcomplex] = ComputeWedgeResponses(audio, 8192);

% figure
% semilogx(fvec, tfmag)

figure
plot(1:impulseLength, audio(1:impulseLength))

[b, a] = deal(zeros(3,numBiquads,numObservations));
for i = 1:numBiquads
    b(:,i,:) = [ones(1,numObservations); -2 * zR(i,:); zR(i,:).^2 + zI(i,:).^2];
    a(:,i,:) = [ones(1,numObservations); -2 * pR(i,:); pR(i,:).^2 + pI(i,:).^2];
end
b(:,1,:) = squeeze(b(:,1,:)) .* k;

input = zeros(3 + sampleDelay, numBiquads + 1);
output = zeros(3, numBiquads);

x = audio;
y = zeros(1, numObservations * step);
for k = 1:numObservations
    for i = 1:step
        idx = (k - 1) * step + i;   % Input sample idx
        input = circshift(input, 1);
        input(1,1) = x(idx);
    
        output = circshift(output, 1);
        for j = 1:numBiquads
            y(idx) = b(1,j,k) * input(sampleDelay + 1,j) + b(2,j,k) * input(sampleDelay + 2,j) + b(3,j,k) * input(sampleDelay + 3,j) - a(2,j,k) * output(2,j) - a(3,j,k) * output(3,j);
            input(sampleDelay + 1,j + 1) = y(idx);
            output(1,j) = y(idx);
        end
    end
end

[ir, tfmag, tfcomplex] = ComputeWedgeResponses(y', 8192);

% figure
% semilogx(fvec, tfmag)

figure
plot(1:numObservations * step, y)

audiowrite("whiteNoise_filt.wav", y, fs)
