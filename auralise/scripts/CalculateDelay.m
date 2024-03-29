%% Calculate the sample delay and fractional delay for a given path length

function [delay, fracDelay, amplitude] = CalculateDelay(pathLength, c, fs, diffPart)
    delay = pathLength * fs / c + 1;     % in seconds // time = distance / speed
    fracDelay = rem(delay, 1);
    delay = floor(delay);
    amplitude = 1 ./ pathLength;

    if nargin > 3
        delayStore = delay;
        delay = delayStore(end, :);
        delay(diffPart) = delayStore(1, diffPart);

        fracDelayStore = fracDelay;
        fracDelay = fracDelayStore(end, :);
        fracDelay(diffPart) = fracDelayStore(1, diffPart);

        amplitudeStore = amplitude;
        amplitude = amplitudeStore(end, :);
        amplitude(diffPart) = amplitudeStore(1, diffPart);
    end
    if pathLength == 0
        amplitude = 1;
    end
end