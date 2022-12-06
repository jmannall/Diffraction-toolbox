%% Plot a spectrogram from a WAV file

function PlotSpectrogramOfWAV(filePath, limits, nfft, savePlot)

    if nargin < 4
        savePlot = false;
    end
    
    if ~exist(filePath, 'file')
        error(['Audio file ', filePath, ' not found'])
    end

    [audio, fs] = audioread(filePath);

    numChannels = size(audio, 2);
    for i = 1:numChannels
        [x, fvec, tvec] = DefaultSpectrogram(audio(:,i), fs, nfft);
    
        fileName = strrep(extractAfter(extractBefore(filePath, '.'), '\'), '_', ' ');
        PlotSpectrogram(x, fvec, tvec, limits, fileName, false, savePlot);
    end
end