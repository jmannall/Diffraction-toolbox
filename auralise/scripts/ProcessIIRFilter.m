function [inputBuffer, outputBuffer, input] = ProcessIIRFilter(input, b, a, inputBuffer, outputBuffer, numBiquads, k)
    inputBuffer = circshift(inputBuffer, 1);
    outputBuffer = circshift(outputBuffer, 1);
    for j = 1:numBiquads
        inputBuffer(1,j) = input;
        input = b(1,j,k) * inputBuffer(1,j) + b(2,j,k) * inputBuffer(2,j) + b(3,j,k) * inputBuffer(3,j) - a(2,j,k) * outputBuffer(2,j) - a(3,j,k) * outputBuffer(3,j);
        outputBuffer(1,j) = input;
    end
end