function [x, fvec, tvec] = DefaultSpectogram(input, nfft, fs)

    windowSize = 400;
    overlap = 300;
    [x, fvec, tvec] = spectrogram(input,windowSize,overlap,nfft,fs,"yaxis");
end