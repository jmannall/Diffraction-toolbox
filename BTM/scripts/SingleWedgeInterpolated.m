function [tfmag, fvec, tfcomplex] = SingleWedgeInterpolated(wedgeLength, wedgeIndex, thetaS, thetaR, radiusS, radiusR, zS, zR, controlparameters, createPlot)
    
    nfft = controlparameters.nfft;
    fs = controlparameters.fs;
    c = controlparameters.c;
    if isfield(controlparameters, 'Rstart')
        controlparameters.Rstart = controlparameters.Rstart - 1 / fs * c;
    end
    input = [1; zeros(5, 1)];
    epsilon = 1e-10;
    if thetaR - thetaS > 180
        % Generate true response
        [~, tfmag, ~, fvec, tfcomplex] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaR, radiusR, radiusS, zS, zR, controlparameters, createPlot);
        
        % Generate reference boundary responses
        [~, tfmagDiffRef, ~, ~, ~] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaS + 180 + epsilon, radiusR, radiusS, zS, zR, controlparameters, createPlot);
        
        zA = CalculateApex(radiusS, radiusR, zS, zR, wedgeLength);
        zA = zA(3);
        pathLength = sqrt(radiusS ^ 2 + (abs(zA - zS) ^ 2)) + sqrt(radiusR ^ 2 + (abs(zA - zR) ^ 2));
        [~, dirIr] = DelayLine(input, pathLength, 3, 1, c, fs);
        %[~, tfmagDirRef, ~, ~, ~] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaS + 180 - epsilon, radiusR, radiusS, zS, zR, controlparameters, createPlot);
        
        % Create scaled response
        DiffRef = tfmagDiffRef.diff1;
        DirRef = IrToTf(dirIr, nfft);
        shift = DirRef - DiffRef;
        scaledResponse = tfmag.diff1 + shift;
        truth = 0;

        % Interpolate between the true and scaled responses
        i = min(1, (thetaR - thetaS - 180) / (wedgeIndex - thetaS - 180 - truth));
        tfmag = (1 - i) * scaledResponse + i * tfmag.diff1;
        
        % Adjust magnitude of tfcomplex
        phase = angle(tfcomplex.diff1);
        tfcomplex = PolarToComplex(10.^(tfmag / 20), phase);
    else
        % Remove diffracted sound if not in the shadow zone 
        [~, dirIr] = DelayLine(input, pathLength, 3, 1, c, fs);
        [tfmag, tfcomplex] = IrToTf(dirIr, nfft);
        fvec = fs/nfft*[0:nfft/2-1];
%         [~, tfmag, ~, fvec, tfcomplex] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaR, radiusR, radiusS, zS, zR, controlparameters, createPlot);
%         tfcomplex = tfcomplex.direct;
%         tfmag = tfmag.direct;
    end
end
