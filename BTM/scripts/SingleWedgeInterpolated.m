function [tfmagOut, fvec, tfcomplex] = SingleWedgeInterpolated(wedgeLength, wedgeIndex, thetaS, thetaR, radiusS, radiusR, zS, zR, controlparameters, createPlot)
    
    nfft = controlparameters.nfft;
    fs = controlparameters.fs;
    c = controlparameters.c;
    L = 0;
    if isfield(controlparameters, 'Rstart')
        controlparameters.Rstart = controlparameters.Rstart - 1 / fs * c;
        L = controlparameters.Rstart;
        if L < 10
            controlparameters = rmfield(controlparameters, 'Rstart');
        end
    end
    epsilon = 1e-10;
    if thetaR - thetaS > 180
        % Removed phase shift for very long paths.
        if L < 10
            input = [1; zeros(5, 1)];
            zA = CalculateApex(radiusS, radiusR, zS, zR, wedgeLength);
            zA = zA(3);

            pathLength = sqrt(radiusS ^ 2 + (abs(zA - zS) ^ 2)) + sqrt(radiusR ^ 2 + (abs(zA - zR) ^ 2));
            [~, dirIr] = DelayLine(input, pathLength, 3, 1, c, fs);
            dirRef = IrToTf(dirIr, nfft);

            scale = pathLength;
        else
            dirRef = zeros(nfft / 2, 1);
            scale = 1;
        end

        % Generate true response
        [~, tfmag, ~, fvec, tfcomplex] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaR, radiusR, radiusS, zS, zR, controlparameters, createPlot);
        
        % Generate reference boundary responses
        [~, tfmagDiffRef, ~, ~, ~] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaS + 180 + epsilon, radiusR, radiusS, zS, zR, controlparameters, createPlot);
        
        %[~, tfmagDirRef, ~, ~, ~] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaS + 180 - epsilon, radiusR, radiusS, zS, zR, controlparameters, createPlot);
        
        % Create scaled response
        diffRef = tfmagDiffRef.diff1;
        shift = dirRef - diffRef;
        scaledResponse = tfmag.diff1 + shift;
        truth = 0;

        % Interpolate between the true and scaled responses
        i = min(1, (thetaR - thetaS - 180) / (wedgeIndex - thetaS - 180 - truth));
        tfmagOut = (1 - i) * scaledResponse + i * tfmag.diff1;
        
        % Adjust magnitude of tfcomplex
        phase = angle(tfcomplex.diff1);
        tfcomplex = scale * PolarToComplex(10.^(tfmagOut / 20), phase);
        tfmagOut = mag2db(abs(tfcomplex));
    else
        % Remove diffracted sound if not in the shadow zone
        l = sqrt(radiusS .^ 2 + zS .^ 2);
        m = sqrt(radiusR .^ 2 + zR .^ 2);

        pathLength = sqrt(l .^ 2 + m .^ 2 - 2 * l .* m .* cosd(thetaR - thetaS));

        input = [1; zeros(11, 1)];
        [~, dirIr] = DelayLine(input, pathLength, 12, 1, c, fs);
        [tfmagOut, tfcomplex] = IrToTf(dirIr, nfft);
        fvec = fs/nfft*[0:nfft/2-1];
%         [~, tfmag, ~, fvec, tfcomplex] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaR, radiusR, radiusS, zS, zR, controlparameters, createPlot);
%         tfcomplex = tfcomplex.direct;
%         tfmag = tfmag.direct;
    end
end
