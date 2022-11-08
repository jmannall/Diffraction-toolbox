function [tfmag, fvec, tfcomplex] = SingleWedgeInterpolated(wedgeLength, wedgeIndex, thetaS, thetaR, radiusS, radiusR, zS, zR, controlparameters, createPlot)
    
    nfft = controlparameters.nfft;
    fs = controlparameters.fs;
    c = controlparameters.c;
    if isfield(controlparameters, 'Rstart')
        controlparameters.Rstart = controlparameters.Rstart - 1 / fs * c;
    end
    epsilon = 1e-10;
    if thetaR - thetaS > 180
        % Generate true response
        [~, tfmag, ~, fvec, tfcomplex] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaR, radiusR, radiusS, zS, zR, controlparameters, createPlot);
        
        % Generate reference boundary responses
        [~, tfmagDiffRef, ~, ~, ~] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaS + 180 + epsilon, radiusR, radiusS, zS, zR, controlparameters, createPlot);
        controlparameters.difforder = 0;
        input = [1; zeros(5, 1)];
        pathLength = sqrt((radiusS + radiusR) ^ 2 + (zR - zS) ^ 2);
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
        [~, tfmag, ~, fvec, tfcomplex] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaR, radiusR, radiusS, zS, zR, controlparameters, createPlot);
        tfcomplex = tfcomplex.direct;
        tfmag = tfmag.direct;
    end
end
