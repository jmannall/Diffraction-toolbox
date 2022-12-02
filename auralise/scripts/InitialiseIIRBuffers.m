%% Initialise buffers for 2nd order biquads in series

function [inputBuffer, outputBuffer] = InitialiseIIRBuffers(numCoeff, numBiquads)

    [inputBuffer, outputBuffer] = deal(zeros(numCoeff, numBiquads));
end