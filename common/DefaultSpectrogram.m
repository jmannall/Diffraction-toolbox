function [x, fvec, tvec] = DefaultSpectrogram(input, fs, nfft)

    if argin < 3
        nfft = 2048;
    end
    windowSize = 400;
    overlap = 300;
    [x, fvec, tvec] = spectrogram(input,windowSize,overlap,nfft,fs,"yaxis");
end