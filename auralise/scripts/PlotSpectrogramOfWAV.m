%% Plot a spectrogram from a WAV file

function PlotSpectrogramOfWAV(filePath, limits, nfft)

    if ~exist(filePath, 'file')
        error(['Audio file ', filePath, ' not found'])
    end

    [audio, fs] = audioread(filePath);

    [x, fvec, tvec] = DefaultSpectrogram(audio, fs, nfft);

    fileName = strrep(extractAfter(extractBefore(filePath, '.'), '\'), '_', ' ');
    PlotSpectrogram(x, fvec, tvec, limits, fileName, false, false);
end