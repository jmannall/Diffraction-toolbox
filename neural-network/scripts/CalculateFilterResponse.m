function [tfmag, fvec, tfcomplex] = CalculateFilterResponse(b, a, nfft, fs)

    x = fft(b, nfft);
    y = fft(a, nfft);
    
    x = squeeze(prod(x,2));
    y = squeeze(prod(y,2));
    
    epsilon = 1e-12;
    F = x ./ (y + epsilon);
    %F = x ./ y;
    tfcomplex = F(1:(end / 2),:);
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