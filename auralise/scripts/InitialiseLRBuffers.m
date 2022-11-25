%% Initialise buffers for a Linkwitz-Riley filterbank

function [inputBuffer, outputBuffer] = InitialiseLRBuffers()

    numBiquads = 4;
    numBands = 4;
    [inputBuffer, outputBuffer] = deal(zeros(3, numBiquads, numBands));
end