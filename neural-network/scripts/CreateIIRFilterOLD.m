function [tfmag, fvec] = CreateIIRFilterOLD(z, p, k, fs)
    
    b = [k, -k * sum(z), k * prod(z)];
    a = [1, -sum(p), prod(p)];
    
    [~, f] = freqz(b, a, 2048);
    b = stripdims(b);
    a = stripdims(a);
    x = fft(b, 4096);
    y = fft(a, 4096);

    tfmag = (20 * log(abs(x ./ y)) / log(10))';
    tfmag = tfmag(1:end / 2);
    
    fvec = (f * fs) / (2 * pi);   
end