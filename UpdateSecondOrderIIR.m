function UpdateSecondOrderIIR(S)

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
    
    y = 1 / (10 ^ 0.3 * S.k ^ 2);
    a = (y * (S.p(1)^2 + 1) - S.z(1)^2 - 1) / (2 * (S.p(1) * y - S.z(1)));
    b = sqrt(1 - a^2);
    z = a+b*1i;
    omegac = angle(z);
    S.fc = omegac / (2 * pi) * S.fs;

    Error = S.k * (Multiply(z, S.z) / Multiply(z, S.p));
    S.Error = 20*log10(Error) + 3;

    hshDCGain = S.k * (Multiply(1, S.z(2)) / Multiply(1, S.p(2)));
    hshNyqGain = S.k * (Multiply(-1, S.z(2)) / Multiply(-1, S.p(2)));

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
    S.dcLabel = uicontrol('Style','text','String',['DC gain: ', num2str(S.dc)],'Position',[1250 760 500 40],'fontsize',20,'HorizontalAlignment','left');
    S.nyqLabel = uicontrol('Style','text','String',['nyq gain: ', num2str(S.nyq)],'Position',[1250 680 500 40],'fontsize',20,'HorizontalAlignment','left');
    S.B0Label = uicontrol('Style','text','String',['B0: ', num2str(S.B0)],'Position',[1250 600 500 40],'fontsize',20,'HorizontalAlignment','left');
    S.fcLabel = uicontrol('Style','text','String',['LPF cut-off frequency: ', num2str(S.fc)],'Position',[1250 520 500 40],'fontsize',20,'HorizontalAlignment','left');
    S.f1Label = uicontrol('Style','text','String',['Hsh f1: ', num2str(S.f1)],'Position',[1250 440 500 40],'fontsize',20,'HorizontalAlignment','left');
    S.f2Label = uicontrol('Style','text','String',['Hsh f2: ', num2str(S.f2)],'Position',[1250 360 500 40],'fontsize',20,'HorizontalAlignment','left');
    S.B1Label = uicontrol('Style','text','String',['Hsh B1: ', num2str(S.B1)],'Position',[1250 280 500 40],'fontsize',20,'HorizontalAlignment','left');

    [y, ~] = S.fcn(S.z, S.p, S.k);   % @(z, p, k) IIRFilter(z,p,k,fs);
    figure(S.fh);
    SetString(S.zLabel,'Zero',S.z);
    SetString(S.pLabel,'Pole',S.p);
    SetString(S.kLabel,'Gain',S.k);
    set(S.LN, 'YData', y);

end