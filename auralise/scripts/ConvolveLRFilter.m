function [output, tfmag, tfcomplex] = ConvolveLRFilter(audio, freqResponse, fs, nfft)

    crossFilt = crossoverFilter( ...
        'NumCrossovers', 3, ...
        'CrossoverFrequencies', [250, 1000, 4000], ...
        'CrossoverSlopes', 24, ...
        'SampleRate', fs);
    [y1,y2,y3,y4] = crossFilt(audio);

    freqResponse = 10 .^ (freqResponse / 20);
    output = freqResponse(1) * y1 + freqResponse(2) * y2 + freqResponse(3) * y3 + freqResponse(4) * y4;

    if nargin < 4
        nfft = 2048;
        disp('No nfft provided. Using default nfft: 2048')
    end
    [tfmag, tfcomplex] = IrToTf(output, nfft);
end