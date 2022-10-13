%% Calculate the sample delay and fractional delay for a given path length

function [delay, fracDelay, amplitude] = CalculateDelay(pathLength, c, fs)
    delay = pathLength * fs / c + 1;     % in seconds // time = distance / speed
    fracDelay = rem(delay, 1);
    delay = floor(delay);
    amplitude = 1 ./ pathLength;
end