function [tf, fvec] = IIRFilter(z, p, k, fs)

[b, a] = zp2tf(z, p, k);

[h, f] = freqz(b, a, 2048);

fvec = (f * fs) / (2 * pi);

tf = 20*log10(abs(h));

end