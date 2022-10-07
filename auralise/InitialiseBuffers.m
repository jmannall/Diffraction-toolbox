function [buffer, read, write, window, overlap, numBuffers, inputBuffer, output] = InitialiseBuffers(delay, samplesPerUpdate, audio, pathLength)

    bufferLength = max(max(delay), samplesPerUpdate);
    buffer = zeros(bufferLength, 1);
    read = 0;
    write = mod(delay(1), bufferLength);
    window = hanning(samplesPerUpdate);
    numBuffers = length(pathLength);
    overlap = floor(samplesPerUpdate / 2);
    inputBuffer = zeros(2, 1);
    output = zeros(length(audio), 1);
end