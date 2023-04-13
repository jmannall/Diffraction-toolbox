%% Create IIR filter diffraction model Kirsch and Ewert
% Currently only valid for zA on the edge
function [tfmag, fvec, tfmagRef, H] = SingleIIRWedge(wedgeLength, wedgeIndex, thetaS, thetaR, rS, rR, zS, zR, n, controlparameters)

    wedgeIndex = deg2rad(wedgeIndex);
    thetaS = deg2rad(thetaS);
    thetaR = deg2rad(thetaR);

    gParameters = struct('wL', wedgeLength, 'wI', wedgeIndex, 'thetaS', thetaS, 'thetaR', thetaR, 'rS', rS, 'rR', rR, 'zS', zS, 'zR', zR);
    if n < 3
        fMin = 20;
    else
        fMin = 10;
    end
    fMax = controlparameters.fs * min(n / 4, 1);

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
    for i = 1:n
        Hi{i} = HshFilter(fvec, fsh(i), gsh(i));
        H = H .* Hi{i};
    end

%     figure
%     semilogx(ft, mag2db(gt), 'x')
%     hold on
%     semilogx(fRef, mag2db(abs(HRef)))
%     semilogx(fi, mag2db(abs(gi)), 'o')
%     semilogx(fsh, mag2db(abs(gi)), 'x')
%     for i = 1:n
%         semilogx(fvec, mag2db(gt(i) * abs(Hi{i})))
%     end
%     grid on
%     xlim([20 22.05e3])
%     ylim([-24 0])
%     semilogx(fvec, mag2db(gt(1) * abs(H)), '--')
%     yticks([-24:6:0])
%     grid on
%     ylim([-24 0])

    dZ = abs(zR - zS);
    pathLength = sqrt((rS + rR) .^ 2 + dZ .^ 2);
    tfmag = mag2db(gt(1) * abs(H) / pathLength);
    tfmagRef = mag2db(abs(HRef) / pathLength);

%     input = [1; zeros(1e4, 1)];
%     [~, dirIr] = DelayLine(input, pathLength, 3, true(size(pathLength)), controlparameters.c, controlparameters.fs);
%     dirRef = IrToTf(dirIr, nfft);
    %tfmag = tfmag + dirRef;
end

