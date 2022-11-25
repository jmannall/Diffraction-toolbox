close all

fs = 48e3;
fc = [250 1000 4000];
[b, a] = CalculateLRCoefficients(fc, fs);

crossFilt = crossoverFilter( ...
    'NumCrossovers', 3, ...
    'CrossoverFrequencies', [250, 1000, 4000], ...
    'CrossoverSlopes', 24, ...
    'SampleRate', fs);

visualize(crossFilt)
[audio, fs] = LoopAudio('audio/whiteNoise96k.wav', 0.5);
audio = [1; zeros(1000, 1)];
fs = 48e3;
crossFilt = crossoverFilter( ...
    'NumCrossovers', 3, ...
    'CrossoverFrequencies', [250, 1000, 4000], ...
    'CrossoverSlopes', 24, ...
    'SampleRate', fs);

% crossFilt = crossoverFilter( ...
%     'NumCrossovers', 2, ...
%     'CrossoverFrequencies', [250, 4000], ...
%     'CrossoverSlopes', 24, ...
%     'SampleRate', fs);

% crossFilt = crossoverFilter( ...
%     'NumCrossovers', 1, ...
%     'CrossoverFrequencies', [250], ...
%     'CrossoverSlopes', 24, ...
%     'SampleRate', fs);


visualize(crossFilt)
cost(crossFilt)

[y1,y2,y3,y4] = crossFilt(audio);

recombine = y1 + 0.5 * y2 + 0.4 * y3 + 0.2 * y4;

nfft = 4096;
fvec = fs/nfft*[0:nfft/2-1];
tfaudio = IrToTf(audio, nfft);
tfrecombine = IrToTf(recombine, nfft);
figure
semilogx(fvec, tfaudio)
hold on
semilogx(fvec, tfrecombine)
legend('audio', 'filters')
xlim([20 20e3])

figure
spectrogram(recombine,120,100,6000,fs,'yaxis')

figure
spectrogram(audio,120,100,6000,fs,'yaxis')

figure('Position',[100,100,800,700])

subplot(5,1,1)
spectrogram(audio,120,100,6000,fs,'yaxis')
title('Noise')

subplot(5,1,2)
spectrogram(y1,120,100,6000,fs,'yaxis')
title('y1')

subplot(5,1,3)
spectrogram(y2,120,100,6000,fs,'yaxis')
title('y2')

subplot(5,1,4)
spectrogram(y3,120,100,6000,fs,'yaxis')
title('y3')

subplot(5,1,5)
spectrogram(y4,120,100,6000,fs,'yaxis')
title('y3')