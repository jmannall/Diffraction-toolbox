%% Increment circular read write index

function [read, write] = IncrementReadWrite(read, write, bufferLength)
        read = mod(read, bufferLength) + 1;
        write = mod(write, bufferLength) + 1;
end