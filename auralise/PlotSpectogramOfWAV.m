function PlotSpectogramOfWAV(filePath, limits, nfft)

    if ~exist(filePath, 'file')
        error(['Audio file ', filePath, ' not found'])
    end

    [audio, fs] = audioread(filePath);

    [x, fvec, tvec] = DefaultSpectogram(audio, nfft, fs);

    fileName = strrep(extractAfter(extractBefore(filePath, '.'), '\'), '_', ' ');
    PlotSpectogram(x, fvec, tvec, limits, fileName, false);
end