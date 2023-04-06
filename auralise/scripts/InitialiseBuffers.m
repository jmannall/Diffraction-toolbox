%% Initialise a buffer (and other data) for a delay line

function [buffer, read, write, window, overlap, numBuffers, inputBuffer, output, windowLength, audio] = InitialiseBuffers(delay, windowLength, audio, pathLength)

    numBuffers = size(pathLength, 2);    % Is delay the same length
    overlap = floor(windowLength / 2);
    inputBuffer = zeros(2, 1);
    outputLength = (numBuffers + 1) * ceil(windowLength / 2) + 1;
    output = zeros(outputLength, 1);
    bufferLength = max(max(delay), windowLength) + overlap;
    buffer = zeros(bufferLength, 1);
    read = 0;
    write = mod(delay(1), bufferLength);
    if numBuffers == 1
        window = [ones(max(windowLength, bufferLength), 1); 0];
        windowLength = length(window) - 1;
        output = zeros(windowLength + 1, 1);
        audio = [audio; zeros(windowLength - length(audio), 1)];
   else
        window = [hanning(windowLength); 0];
        audio = [audio; zeros(outputLength - length(audio), 1)];
    end
end