%% Implement delay line with an iir filter structure

function [output, ir] = DelayLineIIRFilter(audio, pathLength, windowLength, validPath, b, a, c, fs, doFilter)

    b = extractdata(b);
    a = extractdata(a);

    [delay, fracDelay, amplitude] = CalculateDelay(pathLength, c, fs);

    [buffer, read, write, window, overlap, numBuffers, inputBuffer, output] = InitialiseBuffers(delay, windowLength, audio, pathLength);
    maxDelay = max(delay);
    iirLength = 1e3;
    ir = zeros(maxDelay + iirLength + 1, numBuffers);

    [numCoeff, numFilters, ~] = size(b);
    [inputBufferIIR, outputBufferIIR] = InitialiseIIRBuffers(numCoeff, numFilters);     % IIR filter input/output buffers
    
    disp('Process iir filter and delay line')
    for k = 1:numBuffers
        [read, write] = UpdateReadWrite(read, write, overlap, delay, k);
        for i = 1:windowLength
            [idx, read, write, inputBuffer, input] = ProcessSampleData(windowLength, overlap, read, write, inputBuffer, amplitude, validPath, buffer, k, i);
            if doFilter(k)
                [inputBufferIIR, outputBufferIIR, input] = ProcessIIRFilter(input, b(:,:,k), a(:,:,k), inputBufferIIR, outputBufferIIR);
            end
            [inputBuffer, output, buffer] = ProcessSampleOutput(audio, buffer, input, inputBuffer, fracDelay, window, write, idx, k, i, output);
        end
        irIn = zeros(maxDelay + iirLength + 1, 1);
        irIn(delay(k)) = validPath(k) * amplitude(k) * (1 - fracDelay(k));
        irIn(delay(k) + 1) = validPath(k) * amplitude(k) * fracDelay(k);
        if doFilter(k)
            [tempInputBufferIIR, tempOutputBufferIIR] = InitialiseIIRBuffers(numCoeff, numFilters);
            for i = delay(k):delay(k) + iirLength
                [tempInputBufferIIR, tempOutputBufferIIR, ir(i, k)] = ProcessIIRFilter(irIn(i), b(:,:,k), a(:,:,k), tempInputBufferIIR, tempOutputBufferIIR);
            end
        else
            ir(:,k) = irIn;
        end
    end
end