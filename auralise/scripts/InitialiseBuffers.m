%% Initialise a buffer (and other data) for a delay line

function [buffer, read, write, window, overlap, numBuffers, inputBuffer, output] = InitialiseBuffers(delay, windowLength, audio, pathLength)

    bufferLength = max(max(delay), windowLength);
    buffer = zeros(bufferLength, 1);
    read = 0;
    write = mod(delay(1), bufferLength);
    window = hanning(windowLength);
    numBuffers = length(pathLength);    % Is delay the same length
    overlap = floor(windowLength / 2);
    inputBuffer = zeros(2, 1);
    output = zeros(length(audio), 1);
end