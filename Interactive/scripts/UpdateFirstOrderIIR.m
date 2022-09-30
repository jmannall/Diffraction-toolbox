function UpdateFirstOrderIIR(S)
    [y, ~] = S.fcn(S);   % @(S) IIRFilter(S);
    figure(S.fh);
    SetString(S.zLabel, 'Zero: ', S.z);
    SetString(S.pLabel, 'Pole: ', S.p);
    SetString(S.kLabel, 'Gain: ', S.k);
    set(S.LN, 'YData', y);

    dc = abs(S.k * ((1-S.z) / (1-S.p)));
    S.dc = 20*log10(dc);
    nyq = abs(S.k * ((1+S.z) / (1+S.p)));
    S.nyq = 20*log10(nyq);
    B0 = (1 - S.p) * (1 + S.z) / ((1 + S.p) * (1 - S.z));
    S.B0 = 20*log10(B0);
    target = S.dc - 3;
    y = (10^(target / 10) / S.k^2);
    a = (y * (S.p^2 + 1) - S.z^2 - 1) / (2 * (S.p * y - S.z));
    b = sqrt(1 - a^2);
    z = a+b*1i;
    omegac = angle(z);
    S.fcL = omegac / (2 * pi) * S.fs;
    target = S.nyq - 3;
    y = (10^(target / 10) / S.k^2);
    a = (y * (S.p^2 + 1) - S.z^2 - 1) / (2 * (S.p * y - S.z));
    b = sqrt(1 - a^2);
    z = a+b*1i;
    omegac = angle(z);
    S.fcH = omegac / (2 * pi) * S.fs;

    SetString(S.dcLabel, 'DC gain: ', S.dc);
    SetString(S.nyqLabel, 'nyq gain: ', S.nyq);
    SetString(S.shelfLabel, 'Shelf gain: ', S.B0);
    SetString(S.fcLLabel, 'Low cut-off frequency: ', S.fcL);
    SetString(S.fcHLabel, 'High cut-off frequency: ', S.fcH);
end