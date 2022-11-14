function y = ConvolveTfcomplex(x, tfcomplex, numPos)
    % tfcomplex = tfcomplex(:,end) .* ones(size(tfcomplex));
    H = ([tfcomplex; flipud(conj(tfcomplex(2:end,:)))]);
%     h = ifft(H);
%     h = h(1:round(length(H) / 2));
%     % Code to perform Convolution using Overlap Add Method
%     n1 = length(x);
%     n2 = length(h);
%     N = n1+n2-1;
%     y = zeros(1,N);
%     h1 = [h zeros(1,n2-1)];
%     n3 = length(h1);
%     y = zeros(N+n3-n2,1);
%     x = x';
% 
%     window = hanning(n3);
%     overlap = n2 - 1;
% 
%     numPerPos = floor((n1 - n2) / numPos);
%     for i = 1:n2:n1 - n2
%         if i<=(n1-n2-1)
%             x1 = [x(i:i+n3-n2); zeros(n3-n2,1)];
%             %x1 = x(i:i+n3 - 1);
%         else
%             x1 = [x(i:n1); zeros(n3-n2,1)];
%         end
%         x2 = fft(x1);
%         x3 = x2.*H(:,ceil(i / numPerPos));
%         idx = ceil(i / numPerPos);
%         Htest = ifft(H(:,ceil(i / numPerPos)));
%         xTest = conv(x(i:i+n3-n2), Htest);
%         x4 = ifft(x3);
%         if (i==1)
%             y(1:n3) = x4(1:n3);
%         else
%             y(i:i+n3-1) = y(i:i+n3-1) + x4;
%         end
%         figure
%         plot(1:i + n3, y(1:i + n3))
%         hold on
%         %plot(i:i + n3, 0.1 * x(i:i + n3))
%         plot(i:i + n3 - 1, x4)
%         %plot(i:i + length(xTest) - 1, xTest)
%     end


    audioLength = length(x);
    blockSize = 1024;
    numBlocks = ceil(audioLength / blockSize);
    xp = [x'; zeros((numBlocks + 1) * blockSize - audioLength, 1)];

    y = zeros(length(xp),1);
    window = hanning(2 * blockSize);
    for i = 1:numBlocks
        xb = xp((i - 1) * blockSize + 1:(i + 1) * blockSize) .* window;
        u = FftConvolution(xb, H(:,ceil(i / numBlocks * numPos)));
        idx = (i - 1) * blockSize + 1:min((i - 1) * blockSize + length(u), length(xp));
        y(idx) = y(idx) + u(1:length(idx));
%         figure
%         plot(1:idx(end), y(1:idx(end)))
%         hold on
%         plot(idx, u)
    end
end