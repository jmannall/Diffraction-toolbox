%% Process a single pass through a Linkwitz-Riley filter

function [inputBuffer, outputBuffer, output] = ProcessLRFilterBank(input, b, a, inputBuffer, outputBuffer, tfmag, numBands)

    numBiquads = 4;
    freqResponse = 10 .^ (tfmag / 20);

    bandOutputs = zeros(1, numBands);
    for i = 1:numBands
        [inputBuffer(:,:,i), outputBuffer(:,:,i), bandOutputs(i)] = ProcessIIRFilter(input, b(:,:,i), a(:,:,i), inputBuffer(:,:,i), outputBuffer(:,:,i), numBiquads);
    end
    
    output = sum(freqResponse .* bandOutputs);
end