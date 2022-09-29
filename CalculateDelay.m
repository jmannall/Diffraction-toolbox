function [delayLength, fracDelay] = CalculateDelay(pathLength, c, fs)
    delay = pathLength * fs / c;
    delayLength = floor(delay);
    fracDelay = rem(delay, 1);
end