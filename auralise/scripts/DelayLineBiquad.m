%% Implement delay line with a biquad filter structure

function [output, ir] = DelayLineBiquad(audio, pathLength, windowLength, validPath, b, a, numBiquads, c, fs, doBiquad)

    [delay, fracDelay, amplitude] = CalculateDelay(pathLength, c, fs);

    [buffer, read, write, window, overlap, numBuffers, inputBuffer, output] = InitialiseBuffers(delay, windowLength, audio, pathLength);
    maxDelay = max(delay);
    iirLength = 50;
    irIn = zeros(maxDelay + iirLength + 1, 1);
    ir = zeros(maxDelay + iirLength + 1, numBuffers);

    [inputBufferIIR, outputBufferIIR] = InitialiseIIRBuffers(numBiquads);     % IIR filter input/output buffers
    
    disp('Process biquad filter and delay line')
    for k = 1:numBuffers
        [read, write] = UpdateReadWrite(read, write, overlap, delay, k);
        for i = 1:windowLength
            [idx, read, write, inputBuffer, input] = ProcessSampleData(windowLength, overlap, read, write, inputBuffer, amplitude, validPath, buffer, k, i);
            if doBiquad(k)
                [inputBufferIIR, outputBufferIIR, input] = ProcessIIRFilter(input, b(:,:,k), a(:,:,k), inputBufferIIR, outputBufferIIR, numBiquads);
            end
            [inputBuffer, output, buffer] = ProcessSampleOutput(audio, buffer, input, inputBuffer, fracDelay, window, write, idx, k, i, output);
        end
        irIn(delay(k)) = validPath(k) * amplitude(k) * (1 - fracDelay(k));
        irIn(delay(k) + 1) = validPath(k) * amplitude(k) * fracDelay(k);
        if doBiquad(k)
            [tempInputBufferIIR, tempOutputBufferIIR] = InitialiseIIRBuffers(numBiquads);
            for i = delay(k):maxDelay + iirLength
                [tempInputBufferIIR, tempOutputBufferIIR, ir(i, k)] = ProcessIIRFilter(irIn(i), b(:,:,k), a(:,:,k), tempInputBufferIIR, tempOutputBufferIIR, numBiquads);
            end
        else
            ir(:,k) = irIn;
        end

    end
end