%% Normalise audio to peak at 0dB

function [output, fs, savePath] = NormaliseAudio(filePath, nfft)

    % Check audio file exists
    if ~exist(filePath, 'file')
        error(['Audio file ', filePath, ' not found'])
    end

    % Read audio file
    [audio, fs] = audioread(filePath);
    audio = 0.5 * sum(audio, 2);

    % Run spectogram of file
    x = DefaultSpectrogram(audio, fs, nfft);
    x = mag2db(abs(x));

    % Scale audio so the frequency bins peak at 0dB
    gain = 0 - mean(x, 'all');
    gain = db2mag(gain);
    output = gain * audio;

    savePath = [extractBefore(filePath, '.'), '_normalised.wav'];
    audiowrite(savePath, output, fs)
end

