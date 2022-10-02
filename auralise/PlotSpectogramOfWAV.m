function PlotSpectogramOfWAV(fileName, limits, nfft)
    [audio, fs] = audioread(['audio\', fileName]);

    [x, fvec, tvec] = DefaultSpectogram(audio, nfft, fs);
    x = mag2db(abs(x));

    fileName = extractBefore(fileName, ".");
    PlotSpectogram(x, fvec, tvec, limits, fileName);
end