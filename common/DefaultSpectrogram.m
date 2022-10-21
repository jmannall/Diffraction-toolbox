%% Calculate the transfer function, frequency vector and time vector to plot a spectrogram

function [x, fvec, tvec] = DefaultSpectrogram(input, fs, nfft)

    if nargin < 3
        nfft = 2048;
    end
    windowSize = 400;
    overlap = 300;
    [x, fvec, tvec] = spectrogram(input,windowSize,overlap,nfft,fs,"yaxis");
end