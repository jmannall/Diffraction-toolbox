%% Create 1 / n frequency bands from input frequencies 

function [tfmagNBand, fc, fidx] = CreateFrequencyNBands(tfmag, fvec, n)
    df = fvec(2) - fvec(1);
    
    lowFc = 1000; % Start at 1000 Hz and find lowest Fc
    previous = lowFc * 2 ^ (1 / n);
    while (previous-lowFc) > df
        previous = lowFc;
        lowFc = lowFc / 2^(1 / n);
    end
    
    fc = zeros(1, 1000);
    fidx = ones(size(fvec));
    fc(1) = lowFc;
    i = 1;
    while fc(i) <= max(fvec)
        idx = fvec > fc(i);
        fidx(idx) = i;
        i = i + 1;
        fc(i) = 2^(1 / n) * fc(i - 1);
    end
    i = i- 1;
    fc = fc(1:i);
    
    tfmagNBand = CreateNBandMagnitude(tfmag, fidx);
end