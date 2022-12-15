function [b, a] = CalculateLRCoefficients(fc, fs, freqResponse)

    T = 1 / fs;
    K = tan(2 * pi * fc * T / 2);
    K2 = K .^ 2;
    nfft = 8192;
    numCrossovers = length(fc);
    const = ones(1, numCrossovers);
    a1 = 1 + sqrt(2) * K + K2;

    % Lpf
    bLpf(:,1,:) = [K2; 2 * K2; K2] ./ a1;
    bLpf(:,2,:) = bLpf(:,1,:);

    % Hpf
    bHpf(:,1,:) = [const; -2 * const; const] ./ a1;
    bHpf(:,2,:) = bHpf(:,1,:);

    % Lpf and Hpf
    aLpfHpf(:,1,:) = [a1; 2 * (K2 - 1); (1 - K * sqrt(2) + K2)] ./ a1;
    aLpfHpf(:,2,:) = aLpfHpf(:,1,:);

    % Create bands
    bLow= [bLpf(:,:,2), bLpf(:,:,3), bHpf(:,:,3), bLpf(:,:,1)];
    aLow = [aLpfHpf(:,:,2), aLpfHpf(:,:,3), aLpfHpf(:,:,3), aLpfHpf(:,:,1)];

    bMidLow = [bLpf(:,:,2), bLpf(:,:,3), bHpf(:,:,3), bHpf(:,:,1)];

    bMidHigh = [bHpf(:,:,2), bLpf(:,:,1), bHpf(:,:,1), bLpf(:,:,3)];

    bHigh = [bHpf(:,:,2), bLpf(:,:,1), bHpf(:,:,1), bHpf(:,:,3)];
    aHigh = [aLpfHpf(:,:,2), aLpfHpf(:,:,1), aLpfHpf(:,:,1), aLpfHpf(:,:,3)];

    b = cat(3, bLow, bMidLow, bMidHigh, bHigh);
    a = cat(3, aLow, aLow, aHigh, aHigh);

    [tfmag, fvec, ~] = CalculateFilterResponse(b, a, nfft, fs);

    if nargin > 2
        tfmag = tfmag + freqResponse;
    end
    totalTfmag = mag2db(sum(10 .^ (tfmag / 20), 2));

    % Plot bands
    figure
    semilogx(fvec, tfmag)
    hold on
    semilogx(fvec, totalTfmag)
    xlim([20 20e3])
    ylim([-70 10])
    xlabel('Frequency (Hz)')
    ylabel('Magnitude (dB)')
    title('Linkwitz Riley filterbank frequency response')
    grid on
end