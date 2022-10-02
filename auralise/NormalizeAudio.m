function output = NormalizeAudio(filePath, nfft)

    % Read audio file
    [audio, fs] = audioread(filePath);

    % Run spectogram of file
    x = DefaultSpectogram(audio, nfft, fs);
    x = mag2db(abs(x));

    % Scale audio so the frequency bins peak at 0dB
    gain = 0 - max(max(x));
    gain = db2mag(gain);
    output = gain * audio;
end

