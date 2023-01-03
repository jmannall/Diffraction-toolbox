%% Implement delay line with a Linkwitz-Riley filter structure

function [output, ir] = DelayLineLR(audio, pathLength, windowLength, validPath, tfmag, c, fs, doFilter)

    [delay, fracDelay, amplitude] = CalculateDelay(pathLength, c, fs, doFilter);
    amplitude(doFilter) = 1; % Path length 1 / r included in tfmag

    [buffer, read, write, window, overlap, numBuffers, inputBuffer, output, windowLength, audio] = InitialiseBuffers(delay, windowLength, audio, pathLength);
    maxDelay = max(delay);
    iirLength = 1e3;
    ir = zeros(maxDelay + iirLength + 1, numBuffers);

    numBands = 4;
    [inputBufferLR, outputBufferLR] = InitialiseLRBuffers();     % Linwitz-Riley IIR filters input/output buffers
    
    fc = [250 1000 4000];
    
    numFreqResponses = size(tfmag, 1);
    if numFreqResponses == 1
        [b, a] = CalculateLRCoefficients(fc, fs, tfmag);
        tfmag = tfmag .* ones(numBuffers, numBands);
    else
        [b, a] = CalculateLRCoefficients(fc, fs, tfmag(1,:));
    end

    disp('Process L-R filterbank and delay line')
    for k = 1:numBuffers
        [read, write] = UpdateReadWrite(read, write, overlap, delay, k);
        for i = 1:windowLength
            [idx, read, write, inputBuffer, input] = ProcessSampleData(windowLength, overlap, read, write, inputBuffer, amplitude, validPath, buffer, k, i);
            if doFilter(k)
                [inputBufferLR, outputBufferLR, input] = ProcessLRFilterBank(input, b, a, inputBufferLR, outputBufferLR, tfmag(k,:), numBands);
            end
            [inputBuffer, output, buffer] = ProcessSampleOutput(audio, buffer, input, inputBuffer, fracDelay, window, write, idx, k, i, output);
        end
        irIn = zeros(maxDelay + iirLength + 1, 1);
        irIn(delay(k)) = validPath(k) * amplitude(k) * (1 - fracDelay(k));
        irIn(delay(k) + 1) = validPath(k) * amplitude(k) * fracDelay(k);
        if doFilter(k)
            [tempInputBufferLR, tempOutputBufferLR] = InitialiseLRBuffers();
            for i = delay(k):delay(k) + iirLength
                [tempInputBufferLR, tempOutputBufferLR, ir(i, k)] = ProcessLRFilterBank(irIn(i), b, a, tempInputBufferLR, tempOutputBufferLR, tfmag(k,:), numBands);
            end
        else
            ir(:,k) = irIn;
        end
%         bTest = b(:,[1:2,7:8],:);
%         aTest = a(:,[1:2,7:8],:);
%         [test2, fvec] = CalculateFilterResponse(bTest, aTest, 8192, fs);
%         test2 = tfmag(k,:) + test2;
%         test3 = mag2db(sum(10 .^ (test2 / 20), 2));
%         test = IrToTf(ir(:,k), 8192);
%         figure
%         semilogx(fvec, test2)
%         hold on
%         semilogx(fvec, test3)
%         semilogx(fvec, test, '--')
%         xlim([20 20e3])
%         ylim([-70 0])
    end
end