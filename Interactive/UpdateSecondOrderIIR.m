function UpdateSecondOrderIIR(S)

    [y, ~] = S.fcn(S);   % @(S) IIRFilter(S);
    figure(S.fh);
    SetString(S.zLabel,'Zero',S.z);
    SetString(S.pLabel,'Pole',S.p);
    SetString(S.kLabel,'Gain',S.k);
    set(S.LN, 'YData', y);

    dc = 1;
    S.k = dc * Multiply(1, S.p) / Multiply(1, S.z);
%     dc = S.k * (Multiply(1, S.z) / Multiply(1, S.p));

    S.dc = 20*log10(dc);
    nyq = S.k * (Multiply(-1, S.z) / Multiply(-1, S.p));
    S.nyq = 20*log10(nyq);
    low = (2 * pi * 20 / S.fs);
    high = (2 * pi * 20000 / S.fs);
    lB0 = PolarToComplex(1, low);
    hB0 = PolarToComplex(1, high);
    B0 = Multiply(hB0, S.p) * Multiply(lB0, S.z) / (Multiply(lB0, S.p) * Multiply(hB0, S.z));
    S.B0 = 20*log10(B0);
    
    ktemp = dc * Multiply(1, S.p(1)) / Multiply(1, S.z(1));
    y = 1 / (10 ^ 0.3 * ktemp ^ 2);
    a = (y * (S.p(1)^2 + 1) - S.z(1)^2 - 1) / (2 * (S.p(1) * y - S.z(1)));
    b = sqrt(1 - a^2);
    z = a+b*1i;
    omegac = angle(z);
    S.fc = omegac / (2 * pi) * S.fs;

    Error = S.k * (Multiply(z, S.z) / Multiply(z, S.p));
    S.Error = 20*log10(Error) + 3;

    hshDCGain = S.k * (Multiply(1, S.z(2)) / Multiply(1, S.p(2)));
    hshNyqGain = S.k * (Multiply(-1, S.z(2)) / Multiply(-1, S.p(2)));
    
    omegaHigh = 20000 * 2 * pi / S.fs;
    z = cos(omegaHigh) + sin(omegaHigh)*1i;
    highGain = S.k * (Multiply(z, S.z) / Multiply(z, S.p));
    S.highGain = 20*log10(highGain);

    %f1
    y = hshDCGain ^ 2 * 10 ^ 0.3 / S.k ^ 2;
    a = (y * (S.p(2)^2 + 1) - S.z(2)^2 - 1) / (2 * (S.p(2) * y - S.z(2)));
    b = sqrt(1 - a^2);
    z = a+b*1i;
    omegac = angle(z);
    S.f1 = omegac / (2 * pi) * S.fs;
    
    % f2
    y = hshNyqGain ^ 2 / (10 ^ 0.3 * S.k ^ 2);
    a = (y * (S.p(2)^2 + 1) - S.z(2)^2 - 1) / (2 * (S.p(2) * y - S.z(2)));
    b = sqrt(1 - a^2);
    z = a+b*1i;
    omegac = angle(z);
    S.f2 = omegac / (2 * pi) * S.fs;

    % B1
    B1 = Multiply(1, S.p(2)) * Multiply(-1, S.z(2)) / (Multiply(-1, S.p(2)) * Multiply(1, S.z(2)));
    S.B1 = 20*log10(B1);

%     target = S.dc - 3;
%     y = (10^(target / 10) / S.k^2);
%     a = (y * (S.p^2 + 1) - S.z^2 - 1) / (2 * (S.p * y - S.z));
%     b = sqrt(1 - a^2);
%     z = a+b*1i;
%     omegac = angle(z);
%     S.fcL = omegac / (2 * pi) * 48000;
%     target = S.nyq - 3;
%     y = (10^(target / 10) / S.k^2);
%     a = (y * (S.p^2 + 1) - S.z^2 - 1) / (2 * (S.p * y - S.z));
%     b = sqrt(1 - a^2);
%     z = a+b*1i;
%     omegac = angle(z);
%     S.fcH = omegac / (2 * pi) * 48000;
%     
    SetString(S.dcLabel, 'DC gain: ', S.dc);
    SetString(S.nyqLabel, 'nyq gain: ', S.nyq);
    SetString(S.B0Label, 'B0 gain: ', S.B0);
    SetString(S.fcLabel, 'LPF cut-off frequency: ', S.fc);
    SetString(S.f1Label, 'Hsh f1: ', S.f1);
    SetString(S.f2Label, 'Hsh f2: ', S.f2);
    SetString(S.B1Label, 'B1 gain: ', S.B1);
%     set(S.FCy, 'XData', [10, S.fc]);
%     set(S.FCx, 'XData', [S.fc, S.fc]);
%     set(S.gradSlope, 'XData', [S.fc, 20000], 'YData', [-3, S.highGain]);
end