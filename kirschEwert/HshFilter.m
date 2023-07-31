function H = HshFilter(f, fc, g)

    s = (1j * f / fc)';
    H = (1 + sqrt(g) * s) ./ (1 + s / sqrt(g));

%     tfmagRef = mag2db(abs(H));
% 
%     s = 1j * 2 * pi * f;
%     omegaA = 2 * pi * fc;
% 
%     H = (omegaA + sqrt(g) * s) ./ (omegaA + s / sqrt(g));
%     tfmagTest = mag2db(abs(H));
% 
%     fs = 48e3;
%     nfft = 8192;
%     T = 1 / fs;
%     omegaD = (2 / T) * atan(omegaA * T / 2);
%     fcNew = omegaD / (2 * pi);
%     omega = 2 * pi * fc * T;
%     omega = (2 / T) * atan(omega * T / 2);
%     K = omegaD * T * sqrt(g);
% 
%     b = [K + 2 * g; K - 2 * g];
%     a = [K + 2; K - 2];
% 
%     b = b / a(1);
%     a = a / a(1);

%     G = mag2db(g);
    
    %[b, a] = HighShelfCoefficients(fc, G, fs);

%     tfmag = CalculateFilterResponse(b, a, nfft, fs);
    
%     figure
%     semilogx(f, tfmagRef)
%     hold on
%     semilogx(f, tfmag, '--')
%     semilogx(f, tfmagTest, ':')
% 
%     grid on
end