%% Implement delay line with a Linkwitz-Riley filter structure

function output = DelayLineLR(audio, pathLength, windowLength, validPath, tfmag, c, fs)

    [delay, fracDelay, amplitude] = CalculateDelay(pathLength, c, fs);
    amplitude = ones(size(amplitude));

    [buffer, read, write, window, overlap, numBuffers, inputBuffer, output] = InitialiseBuffers(delay, windowLength, audio, pathLength);
    
    numBands = 4;
    [inputBufferLR, outputBufferLR] = InitialiseLRBuffers();     % Linwitz-Riley IIR filters input/output buffers
    
    fc = [250 1000 4000];
    
    numFreqResponses = size(tfmag, 1);
    if numFreqResponses == 1
        [b, a] = CalculateLRCoefficients(fc, fs, tfmag);
        tfmag = tfmag .* ones(numBuffers, numBands);
    else
        [b, a] = CalculateLRCoefficients(fc, fs);
    end

    disp('Process L-R filterbank and delay line')
    for k = 1:numBuffers
        [read, write] = UpdateReadWrite(read, write, overlap, delay, k);
        for i = 1:windowLength
            [idx, read, write, inputBuffer, input] = ProcessSampleData(windowLength, overlap, read, write, inputBuffer, amplitude, validPath, buffer, k, i);
            [inputBufferLR, outputBufferLR, input] = ProcessLRFilterBank(input, b, a, inputBufferLR, outputBufferLR, tfmag(k,:), numBands);
            [inputBuffer, output, buffer] = ProcessSampleOutput(audio, buffer, input, inputBuffer, fracDelay, window, write, idx, k, i, output);
        end
    end
end