function [idx, read, write, inputBuffer, input] = ProcessSampleData(windowLength, overlap, read, write, inputBuffer, amplitude, validPath, buffer, k, i)
    
    bufferLength = length(buffer);
    idx = (k - 1) * (windowLength - overlap) + i;
    [read, write] = IncrementReadWrite(read, write, bufferLength);
    inputBuffer = circshift(inputBuffer, 1);
    input = validPath(k) * amplitude(k) * buffer(read);
end