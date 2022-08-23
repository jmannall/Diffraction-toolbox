function UpdateFirstOrderIIR(S)
    [y, ~] = S.fcn(S.z, S.p, S.k);   % @(z, p, k) IIRFilter(z,p,k,fs);
    figure(S.fh);
    SetString(S.zLabel,'Zero',S.z);
    SetString(S.pLabel,'Pole',S.p);
    SetString(S.kLabel,'Gain',S.k);
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
    
    S.dcLabel = uicontrol('Style','text','String',['DC gain: ', num2str(S.dc)],'Position',[1250 600 500 40],'fontsize',20,'HorizontalAlignment','left');
    S.dcLabel = uicontrol('Style','text','String',['nyq gain: ', num2str(S.nyq)],'Position',[1250 520 500 40],'fontsize',20,'HorizontalAlignment','left');
    S.dcLabel = uicontrol('Style','text','String',['Shelf gain: ', num2str(S.B0)],'Position',[1250 440 500 40],'fontsize',20,'HorizontalAlignment','left');
    S.fcLabel = uicontrol('Style','text','String',['Low cut-off frequency: ', num2str(S.fcL)],'Position',[1250 360 500 40],'fontsize',20,'HorizontalAlignment','left');
    S.fcLabel = uicontrol('Style','text','String',['High cut-off frequency: ', num2str(S.fcH)],'Position',[1250 280 500 40],'fontsize',20,'HorizontalAlignment','left');
end