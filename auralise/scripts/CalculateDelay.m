%% Calculate the sample delay and fractional delay for a given path length

function [delay, fracDelay, amplitude] = CalculateDelay(pathLength, c, fs, diffPart)
    delay = pathLength * fs / c + 1;     % in seconds // time = distance / speed
    fracDelay = rem(delay, 1);
    delay = floor(delay);
    amplitude = 1 ./ pathLength;

    if nargin > 3
        delayStore = delay;
        delay = delayStore(:, 2);
        delay(diffPart) = delayStore(diffPart, 1);

        fracDelayStore = fracDelay;
        fracDelay = fracDelayStore(:, 2);
        fracDelay(diffPart) = fracDelayStore(diffPart, 1);

        amplitudeStore = amplitude;
        amplitude = amplitudeStore(:, 2);
        amplitude(diffPart) = amplitudeStore(diffPart, 1);
    end
end