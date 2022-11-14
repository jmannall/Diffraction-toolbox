%% Initialise a buffer (and other data) for a delay line

function [buffer, read, write, window, overlap, numBuffers, inputBuffer, output, windowLength, audio] = InitialiseBuffers(delay, windowLength, audio, pathLength)

    numBuffers = length(pathLength);    % Is delay the same length
    overlap = floor(windowLength / 2);
    inputBuffer = zeros(2, 1);
    output = zeros(length(audio), 1);
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
    end
end