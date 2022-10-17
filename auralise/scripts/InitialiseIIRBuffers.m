%% Initialise buffers for a 2nd order biquads in series

function [inputBuffer, outputBuffer] = InitialiseIIRBuffers(numBiquads)
    [inputBuffer, outputBuffer] = deal(zeros(3, numBiquads));
end