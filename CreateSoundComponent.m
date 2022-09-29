function [sampleDelay, fracDelay, bufferLength, read, write, buffer] = CreateSoundComponent(pathLength, c, fs, samplesPerUpdate)
    delay = pathLength / c;     % in seconds // time = distance / speed
    sampleDelay = floor(delay * fs + 1);
    fracDelay = rem(delay * fs + 1, 1);
    bufferLength = max(max(sampleDelay), samplesPerUpdate);
    read = 0;
    write = mod(sampleDelay(1), bufferLength);
    buffer = zeros(bufferLength, 1);
end