%% Calculate magnitude response from impulse response

function [tfmag, tfcomplex, ir] = IrToTf(ir, nfft)  
    F = fft(ir,nfft);
    tfcomplex = F(1:nfft/2,:);
    tfmag = 20*log10(abs(F(1:nfft/2,:)));
end