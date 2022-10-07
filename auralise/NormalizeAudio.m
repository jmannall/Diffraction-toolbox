function [output, fs] = NormalizeAudio(filePath, nfft)

    % Check audio file exists
    if ~exist(filePath, 'file')
        error(['Audio file ', filePath, ' not found'])
    end

    % Read audio file
    [audio, fs] = audioread(filePath);

    % Run spectogram of file
    x = DefaultSpectogram(audio, nfft, fs);
    x = mag2db(abs(x));

    % Scale audio so the frequency bins peak at 0dB
    gain = 0 - mean(x, 'all');
    gain = db2mag(gain);
    output = gain * audio;

    savePath = [extractBefore(filePath, '.'), '_Normalised.wav'];
    audiowrite(savePath, output, fs)
end

