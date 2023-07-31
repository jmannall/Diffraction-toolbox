%% Create IIR filter diffraction model Kirsch and Ewert
% Currently only valid for zA on the edge
%function [tfmag, fvec, tfmagRef, H] = SingleIIRWedge(wedgeLength, wedgeIndex, thetaS, thetaR, rS, rR, zS, zR, n, controlparameters)

function [tfmag, fvec] = SingleIIRWedge(wedgeLength, wedgeIndex, thetaS, thetaR, rS, rR, zS, zR, n, controlparameters)

%     fs = controlparameters.fs;
%     nfft = controlparameters.nfft;
%     c = controlparameters.c;
%     T = 1 / fs;
%     
%     [zA, phii] = CalculateApex(rS, rR, zS, zR, 100, true);
%     zA = zA(:,3);
%     
%     bendingAngle = abs(thetaR - thetaS);
%     minAngle = min(thetaS, wedgeIndex - thetaR);
%     bendingAngle = max(0.001, deg2rad(bendingAngle) - pi);  % Adjust bending angle
%     minAngle = deg2rad(minAngle);
%     wedgeIndex = deg2rad(wedgeIndex);
%     phii = deg2rad(phii);
%     
%     dS = sqrt(rS ^ 2 + (zA - zS) ^ 2);
%     dR = sqrt(rR ^ 2 + (zR - zA) ^ 2);
%     d = (2 * dS * dR) / (dS + dR); % 4
%     maxBendingAngle = wedgeIndex - (minAngle + pi);
%     mp = 1 - 0.75 * tanh(1 / (2 * bendingAngle)) * sqrt(tanh(2 * minAngle)); % 9
%     mw = (1 - 0.75 * (bendingAngle / maxBendingAngle) * sqrt(sin(-wedgeIndex / 2)))^(-1); % sin(-wedgeIndex / 2) in the paper. But this leads to complex results. This gives a closer match to UTD too so probably a typo?. -Maybe not varies by case
%     % 9
%     fc = (c * real(mw) * mp) / (3 * pi * d * (1 - cos(bendingAngle)) * sin(phii) ^ 2); % 11
% 
%     %% LPF
%     f0 = 1.11 * fc * (15.6 / fc) ^ 0.141; % 4
%     [bLpf, aLpf] = EDlpf(f0, T); % 5
% 
%     pathLength = dS + dR;
%     k = aLpf(1);
%     aLpf = aLpf / k;
%     bLpf = bLpf / (k * pathLength);
%     
%     %% Hsh
%     actual = (log2(fs / 2) - log2(min(f0, fs / 2))) * 6; % 4
%     target = (log2(fs / 2) - log2(min(fc, fs / 2))) * 3; % 4
%     
%     fsh = 209 * fc * (15.6 / fc) ^ 0.827; % 4
%     G = max(0, actual - target); % 1
%     [bHsh, aHsh] = EDhshOld(fsh, G, T); % 10
%     
%     k = aHsh(1);
%     aHsh = aHsh / k;
%     bHsh = bHsh / k;
% 
%     %% Combine filters
% 
%     b = [bLpf bHsh];
%     a = [aLpf aHsh];
%    [tfmag, fvec] = CalculateFilterResponse(b, a, nfft, fs);

    wedgeIndex = deg2rad(wedgeIndex);
    thetaS = deg2rad(thetaS);
    thetaR = deg2rad(thetaR);

    gParameters = struct('wL', wedgeLength, 'wI', wedgeIndex, 'thetaS', thetaS, 'thetaR', thetaR, 'rS', rS, 'rR', rR, 'zS', zS, 'zR', zR, ...
        'zA', 0, 'phii', 0, 'dS', 0, 'dR', 0, 'd', 0, 't0', 0, 'v', 0);
    if n < 3
        fMin = 20;
    else
        fMin = 10;
    end
    fMax = min(controlparameters.fs, 48e3) * min(n / 4, 1);

    ft = logspace(log10(fMin), log10(fMax), n + 1); % Precompute
    gt = abs(CalculateUDFATarget(ft, gParameters, controlparameters)); % 1 + n -> 337

    fRef = controlparameters.fs/controlparameters.nfft*[0:controlparameters.nfft/2-1];
    HRef = CalculateUDFATarget(fRef, gParameters, controlparameters);

    [fsh, gsh, fi, gi] = deal(zeros(1, n));
    for i = 1:n
        fi(i) = ft(i) * sqrt(ft(i + 1) / ft(i));    % Precompute
        gsh(i) = gt(i + 1) / gt(i); % 1
        gi(i) = abs(CalculateUDFATarget(fi(i), gParameters, controlparameters));
        gd = gi(i) / gt(i); % 1
        fsh(i) = fi(i) * sqrt((gd ^ 2 - gsh(i) ^ 2) / (gsh(i) * (1 - gd ^ 2))) * (1 + (gsh(i) ^ 2) / 12); % 13 -> Think this fixes a typo in eq (20)?
    end

    fs = controlparameters.fs;
    nfft = controlparameters.nfft;
    fvec = fs/nfft*[0:nfft/2-1];
    H = ones(nfft / 2, 1);
    T = 1 / fs;
    Hi = cell(1,n);
    test = cell(1,n);
    b = cell(1,n);
    a = cell(1,n);

    for i = 1:n
        Hi{i} = HshFilter(fvec, fsh(i), gsh(i));
        [b{i}, a{i}] = EDhsh(fsh(i), gsh(i), T);
    end

    for i = 1:n
        test{i} = mag2db(abs(freqz(b{i}, a{i}, fvec, fs)));
        H = H .* Hi{i};
    end

    figure
    semilogx(fvec, test{1})
    hold on
    grid on
    semilogx(fvec, test{2})
    semilogx(fvec, test{3})
    semilogx(fvec, test{4})
    semilogx(fvec, mag2db(gt(1)) + test{1} + test{2} + test{3} + test{4})
    semilogx(fvec, mag2db(gt(1) * abs(H)), '--')
    xlim([20 22.05e3])
    yticks([-24:6:0])
    ylim([-24 0])

    figure
    semilogx(ft, mag2db(gt), 'x')
    hold on
    semilogx(fRef, mag2db(abs(HRef)), '-.')
    semilogx(fi, mag2db(abs(gi)), 'o')
    semilogx(fsh, mag2db(abs(gi)), 'x')
    for i = 1:n
        semilogx(fvec, mag2db(gt(i) * abs(Hi{i})))
        semilogx(fvec, mag2db(gt(i)) + test{i}, '--')
    end
    grid on
    xlim([20 22.05e3])
    ylim([-24 0])
    semilogx(fvec, mag2db(gt(1) * abs(H)))
    semilogx(fvec, mag2db(gt(1)) + test{1} + test{2} + test{3} + test{4}, '--')
    yticks([-24:6:0])
    grid on
    ylim([-24 0])

    dZ = abs(zR - zS);
    pathLength = sqrt((rS + rR) .^ 2 + dZ .^ 2);
    tfmag = mag2db(gt(1) * abs(H) / pathLength);
    %tfmag = mag2db(gt(1) * abs(H));
end

