%% Process a single pass through an IIR filter

function [inputBuffer, outputBuffer, input] = ProcessIIRFilter(input, b, a, inputBuffer, outputBuffer, numBiquads)
    inputBuffer = circshift(inputBuffer, 1);
    outputBuffer = circshift(outputBuffer, 1);
    for j = 1:numBiquads
        inputBuffer(1,j) = input;
        % input = b(:,j) .* inputBuffer(:,j);     % Should be the same as
        % below and epand to any order?
        input = b(1,j) * inputBuffer(1,j) + b(2,j) * inputBuffer(2,j) + b(3,j) * inputBuffer(3,j) - a(2,j) * outputBuffer(2,j) - a(3,j) * outputBuffer(3,j);
        outputBuffer(1,j) = input;
    end
end