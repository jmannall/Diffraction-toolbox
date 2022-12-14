function [output, tfmag, tfcomplex] = DelayLineLRFilter(audio, tfmag, pathLength, windowLength, validPath, c, fs, nfft)
    
%     if length(audio) < windowLength
%         audio = [audio; zeros(windowLength * 100, 1)];
%     end
    
    % [audio, ~] = DelayLine(audio, pathLength, windowLength, validPath, c, fs);

    crossFilt = crossoverFilter( ...
        'NumCrossovers', 3, ...
        'CrossoverFrequencies', [250, 1000, 4000], ...
        'CrossoverSlopes', 24, ...
        'SampleRate', fs);

    [delay, fracDelay, amplitude] = CalculateDelay(pathLength, c, fs);
    const = ones(size(amplitude));

    [buffer, read, write, window, overlap, numBuffers, inputBuffer, output, windowLength, audio] = InitialiseBuffers(delay, windowLength, audio, pathLength);
    ir = zeros(length(buffer), numBuffers);
    
    maxIdx = length(output);
    % disp('Process delay line')
    for k = 1:numBuffers
        ir(delay(k), k) = validPath(k) * amplitude(k) * (1 - fracDelay(k));
        ir(delay(k) + 1, k) = validPath(k) * amplitude(k) * fracDelay(k);
        [read, write] = UpdateReadWrite(read, write, overlap, delay, k);
        for i = 1:windowLength
            [idx, read, write, inputBuffer, input] = ProcessSampleData(windowLength, overlap, read, write, inputBuffer, const, validPath, buffer, k, i);
            [inputBuffer, output, buffer] = ProcessSampleOutput(audio, buffer, input, inputBuffer, fracDelay, window, write, idx, k, i, output);
        end

        idx = (k - 1) * windowLength / 2 + 1:min(k * windowLength, maxIdx);
        [y1,y2,y3,y4] = crossFilt(output(idx));
        freqResponse = 10 .^ (tfmag(k,:) / 20);
        output(idx) = output(idx) + freqResponse(1) * y1 + freqResponse(2) * y2 + freqResponse(3) * y3 + freqResponse(4) * y4;
    end

    if nargin < 4
        nfft = 2048;
        disp('No nfft provided. Using default nfft: 2048')
    end

    [tfmag, tfcomplex] = IrToTf(output, nfft);
end