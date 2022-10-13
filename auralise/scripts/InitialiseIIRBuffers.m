function [inputBuffer, outputBuffer] = InitialiseIIRBuffers(numBiquads)
    [inputBuffer, outputBuffer] = deal(zeros(3, numBiquads));
end