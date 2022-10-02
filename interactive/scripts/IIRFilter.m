function [tfmag, fvec] = IIRFilter(S)

    % Inputs
    z = S.z;
    p = S.p;
    k = S.k;
    fs = S.fs;

    % Create IIR filter and find magnitude response and frequency vector
    [b, a] = zp2tf(z, p, k);
    [h, f] = freqz(b, a, 2048);
    
    tfmag = 20*log10(abs(h));
    fvec = (f * fs) / (2 * pi);
end