function [a_lm, targets, DI, invDF] = CalculateSphericalHarmonics(data, freqBands, fs, normalise)
    
    [N, numPositions] = size(data.IR);
    numFreqBands = length(freqBands);  % Number of frequency bands to analyze
    
    [PHI, THETA] = meshgrid(unique(data.Phi), unique(data.Theta));
    
    PHI = [PHI, PHI(:,1)];
    THETA = [THETA, THETA(:,1)];
    
    phi = PHI(1,:);
    theta = THETA(:,1);
    
    PHI = deg2rad(PHI);
    THETA = deg2rad(THETA);
    
    idx = zeros(size(PHI));
    for i = 1:numPositions
        p = data.Phi(i);
        t = data.Theta(i);
    
        idxP = phi == p;
        idxT = theta == t;
        idx(idxT, idxP) = i;
    end
    idx(1,:) = 1;
    idx(end,:) = numPositions;
    
    X = cos(THETA);
    Y = sin(PHI) .* sin(THETA);
    Z = cos(PHI) .* sin(THETA);
    
    targets = cell(numFreqBands, 1);
    a_lm = cell(numFreqBands, 1);
    
    % Precompute FFT for all impulse responses
    fftResponses = fft(data.IR, N, 1);  % Perform FFT along time axis
    frequencies = (0:N-1) * (fs / N);  % Frequency vector
    
    magnitudeResponse = zeros(numPositions, numFreqBands);
    
    % Loop over frequency bands
    for k = 1:numFreqBands
        % Get the center frequency for the current band
        centerFreq = freqBands(k);
        
        % Define a frequency band around the center frequency (e.g., octave band)
        fLower = centerFreq / sqrt(2);  % Lower bound of the band
        fUpper = centerFreq * sqrt(2);  % Upper bound of the band
        
        % Find indices of frequencies within this band
        bandIdx = frequencies >= fLower & frequencies <= fUpper;
        
        % Compute the magnitude of the frequency response in this band
        magnitudeResponse(:,k) = mean(abs(fftResponses(bandIdx, :)), 1);  % Mean across band
    end
    
    % Normalise
    if (normalise)
        %check = max(sum(magnitudeResponse, 2), [], 'all') * numFreqBands;
        magnitudeResponse = magnitudeResponse / max(sum(magnitudeResponse, 2), [], 'all') * numFreqBands;
    end

    for k = 1:numFreqBands
        store = magnitudeResponse(:,k);
        targets{k} = store(idx);

        %di(k) = 10 * log10(max(targets{k} .^ 2) / mean(targets{k} .^ 2));
        %df(k) = 1 / (10 ^ (di(k) / 20));
    end
    

    %%
    
    % Parameters
    lMax = 6;  % Maximum degree of spherical harmonics (adjust as needed)
    [numTheta, numPhi] = size(THETA);
    %numPhi = 361;      % 361
    filename = [DataHash({lMax, THETA, PHI}) '.mat'];

    if (exist(filename, "file"))
        load(filename, 'Y_lmStore')
    else
        for l = 0:lMax
            for m = -l:l
                % Initialize the sum for this (l, m) pair
                Y_lmStore{l^2 + l + m + 1} = zeros(numTheta, numPhi);
                
                % Calculate spherical harmonic for each grid point
                for i = 1:numTheta
                    for j = 1:numPhi
                        Y_lmStore{l^2 + l + m + 1}(i, j) = harmonicY(l, m, THETA(i, j), PHI(i, j));  % Compute Y_lm(theta, phi)
                    end
                end
            end
        end
        save(filename, 'Y_lmStore', '-mat')
    end    
    
    %%
    % close all
    
    f1 = figure;
    f2 = figure;
    sinTheta = sin(THETA);  % same size as Theta
    for k = 1:numFreqBands
        f = targets{k};
        fReconstructed = zeros(numTheta, numPhi);
        % Loop over all (l, m) pairs to compute spherical harmonic coefficients
        for l = 0:lMax
            for m = -l:l
                % Initialize the sum for this (l, m) pair
                Y_lm = Y_lmStore{l^2 + l + m + 1};
        
                % Compute the integral using numerical integration
                integrand = f .* conj(Y_lm) .* sinTheta;
                % Integrate over theta (rows) and phi (columns) using trapz
                integral_theta = THETA(2,2) * trapz(integrand, 1);  % integrate over theta (dim=1)
                a_lm{k}(l^2 + l + m + 1) = PHI(2,2) * trapz(integral_theta, 2);  % integrate over phi (dim=2)
        
                fReconstructed = fReconstructed + (a_lm{k}(l^2 + l + m + 1) .* Y_lm);
            end
            error = mean(abs(mag2db(f) - mag2db(abs(fReconstructed))), 'all');
            if (error < 0.5)
                disp(['End early: ', num2str(l)]);
                break;
            end
        end
        figure(f1)
        subplot(3, 3, k);
        surf(X, Y, Z, mag2db(abs(fReconstructed)), 'EdgeColor', 'none');  % 3D surface plot
        colormap(jet);  % Color map for intensity
        colorbar;  % Add color bar to indicate magnitude levels
        axis equal;  % Ensure the plot is spherical
        title(sprintf('Directivity at %d Hz', freqBands(k)));
        xlabel('X'); ylabel('Y'); zlabel('Z');
        clim([-40 40])
        view([60 15])
    
        figure(f2)
        subplot(3, 3, k);
        surf(X, Y, Z, mag2db(abs(f)), 'EdgeColor', 'none');  % 3D surface plot
        colormap(jet);  % Color map for intensity
        colorbar;  % Add color bar to indicate magnitude levels
        axis equal;  % Ensure the plot is spherical
        title(sprintf('Directivity at %d Hz', freqBands(k)));
        xlabel('X'); ylabel('Y'); zlabel('Z');
        clim([-40 40])
        view([60 15])

        power = abs(f).^2;        
        DI(k) = 10 * log10( max(power, [], 'all') / mean(power, 'all'));
        invDF(k) = 1 / db2mag(DI(k));
    end
end

