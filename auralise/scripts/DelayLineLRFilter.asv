function [output, tfmag, tfcomplex] = DelayLineLRFilter(audio, freqResponse, pathLength, windowLength, validPath, c, fs, nfft)
    
%     if length(audio) < windowLength
%         audio = [audio; zeros(windowLength * 100, 1)];
%     end
    
    [audio, ~] = DelayLine(audio, pathLength, windowLength, validPath, c, fs);

    crossFilt = crossoverFilter( ...
        'NumCrossovers', 3, ...
        'CrossoverFrequencies', [250, 1000, 4000], ...
        'CrossoverSlopes', 24, ...
        'SampleRate', fs);
    [y1,y2,y3,y4] = crossFilt(audio);

    freqResponse = pathLength * 10 .^ (freqResponse / 20);
    output = freqResponse(1) * y1 + freqResponse(2) * y2 + freqResponse(3) * y3 + freqResponse(4) * y4;
    if nargin < 4
        nfft = 2048;
        disp('No nfft provided. Using default nfft: 2048')
    end

    [tfmag, tfcomplex] = IrToTf(output, nfft);
end