%% Implement delay line

function [output, ir] = DelayLine(audio, pathLength, validPath, c, samplesPerUpdate, fs)

    [delay, fracDelay, amplitude] = CalculateDelay(pathLength, c, fs);

    [buffer, read, write, window, overlap, numBuffers, inputBuffer, output] = InitialiseBuffers(delay, samplesPerUpdate, audio, pathLength);
    ir = zeros(length(buffer), numBuffers);

    disp('Process delay line')
    for k = 1:numBuffers
        ir(delay(k), k) = validPath(k) * amplitude(k) * (1 - fracDelay(k));
        ir(delay(k) + 1, k) = validPath(k) * amplitude(k) * fracDelay(k);
        [read, write] = UpdateReadWrite(read, write, overlap, delay, k);
        for i = 1:samplesPerUpdate
            [idx, read, write, inputBuffer, input] = ProcessSampleData(samplesPerUpdate, overlap, read, write, inputBuffer, amplitude, validPath, buffer, k, i);
            [inputBuffer, output, buffer] = ProcessSampleOutput(audio, buffer, input, inputBuffer, fracDelay, window, write, idx, k, i, output);
        end
    end
end