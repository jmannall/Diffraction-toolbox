function output = ConvolveIR(audio, ir, validPath, samplesPerUpdate)

    bufferLength = size(ir, 1);
    buffer = zeros(bufferLength, 1);
    write = 0;
    read = bufferLength - 1:-1:0;
    window = hanning(samplesPerUpdate);
    overlap = floor(samplesPerUpdate / 2);
    numBuffers = size(ir, 2);
    output = zeros(length(audio), 1);
    
    disp('Convolve impulse response')
    for k = 1:numBuffers
        read = read - overlap;
        write = write - overlap;
        for i = 1:samplesPerUpdate
            idx = (k - 1) * (samplesPerUpdate - overlap) + i;
            [read, write] = IncrementReadWrite(read, write, bufferLength);
            input = validPath(k) * buffer(read)' * ir(:, k);
            output(idx) = output(idx) + window(i) * input;
            buffer(write) = audio(idx);
        end
    end
end