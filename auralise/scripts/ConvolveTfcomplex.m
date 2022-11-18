function output = ConvolveTfcomplex(audio, tfcomplex, numPos)
    % tfcomplex = tfcomplex(:,end) .* ones(size(tfcomplex));
    H = ([tfcomplex; flipud(conj(tfcomplex(2:end,:)))]);

    audioLength = length(audio);
    blockSize = 1024;
    numBlocks = ceil(audioLength / blockSize);
    xp = [audio; zeros((numBlocks + 1) * blockSize - audioLength, 1)];

    output = zeros(length(xp),1);
    window = hanning(2 * blockSize);
    for i = 1:numBlocks
        xb = xp((i - 1) * blockSize + 1:(i + 1) * blockSize) .* window;
        u = FftConvolution(xb, H(:,ceil(i / numBlocks * numPos)));
        idx = (i - 1) * blockSize + 1:min((i - 1) * blockSize + length(u), length(xp));
        output(idx) = output(idx) + u(1:length(idx));
    end
end