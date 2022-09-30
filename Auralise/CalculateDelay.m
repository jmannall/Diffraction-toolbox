%% Calculate the sample delay and fractional delay for a given path length

function [sampleDelay, fracDelay] = CalculateDelay(pathLength, c, fs)
    sampleDelay = pathLength * fs / c;
    fracDelay = rem(sampleDelay, 1);
    sampleDelay = floor(sampleDelay);
end