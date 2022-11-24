%% Calculate the transfer function and magnitude response of an impulse response

function [tfmag, tfcomplex, ir] = IrToTf(ir, nfft)  
    F = fft(ir,nfft);
    tfcomplex = F(1:nfft/2,:);
    tfmag = 20*log10(max(abs(F(1:nfft/2,:)), 1e-10));
end