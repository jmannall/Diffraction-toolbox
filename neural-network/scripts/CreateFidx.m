function fidx = CreateFidx(fvec, n)    
    % idxLo = find(fvec > 20, 1);
    % idxHi = find(fvec < 20e3, 1, "last");
    idxLo = 1;
    idxHi = length(fvec);

    fvec = fvec(idxLo:idxHi);

    df = fvec(2) - fvec(1);
    
    lowFc = 1000; % Start at 1000 Hz and find lowest Fc
    previous = lowFc * 2 ^ (1 / n);
    while (previous-lowFc) > df
        previous = lowFc;
        lowFc = lowFc / 2^(1 / n);
    end
    
    fc = zeros(1, 1000);
    fc(1) = lowFc;
    i = 1;
    while fc(i) <= max(fvec)
        i = i + 1;
        fc(i) = 2^(1 / n) * fc(i - 1);
    end
    numFreq = i;
    fc = fc(1:numFreq);
    fLow = fc / sqrt(2 ^ (1 / n));
    
    fidx = ones(size(fvec));
    for i = 2:numFreq
        idx = fvec > fLow(i);
        fidx(idx) = i;
    end
end