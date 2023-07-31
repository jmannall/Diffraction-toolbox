function [tfmag, fvec, tfcomplex] = CalculateFilterResponseLR(b, a, nfft, fs)

    x = fft(b, nfft);
    y = fft(a, nfft);
     
    for i = 1:4
        xStore = squeeze(x(:,2 * i - 1,:) .* x(:,2 * i,:));
        yStore = squeeze(y(:,2 * i - 1,:) .* y(:,2 * i,:));

        F = xStore ./ yStore;
        tf(:,i,:) = F(1:(end / 2),:);
    end

    epsilon = 1e-12;
    tfcomplex = squeeze(tf(:,1,:) .* (tf(:,2,:) + tf(:,3,:)) .* tf(:,4,:));
    tfmag = (20*log(abs(tfcomplex)+epsilon) / log(10));

    %% Create fvec 

    if nargin == 6
        fvec = 0;
        disp('No sample rate provided for fvec');
    else
        [~, f] = freqz(b(:,1), a(:,1), nfft / 2);
        fvec = (f * fs) / (2 * pi);
    end
end