%% Convolve impulse responses with input audio

function output = ConvolveIR(audio, ir, windowLength, validPath)

    bufferLength = size(ir, 1);
    buffer = zeros(bufferLength, 1);
    write = 0;
    read = bufferLength - 1:-1:0;
    window = hanning(windowLength);
    overlap = floor(windowLength / 2);
    numBuffers = size(ir, 2);
    output = zeros(length(audio), 1);
    
    % Default validPath
    if nargin < 4
        validPath  = ones(1, numBuffers);
    end

    disp('Convolve impulse response')
    for k = 1:numBuffers
        read = read - overlap;
        write = write - overlap;
        for i = 1:windowLength
            idx = (k - 1) * (windowLength - overlap) + i;
            [read, write] = IncrementReadWrite(read, write, bufferLength);
            input = validPath(k) * buffer(read)' * ir(:, k);
            output(idx) = output(idx) + window(i) * input;
            buffer(write) = audio(idx);
        end
    end
end