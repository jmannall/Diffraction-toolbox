function [tfmag, fvec, tfcomplex] = SingleWedgeInterpolated(wedgeLength, wedgeIndex, thetaS, thetaR, radiusS, radiusR, zS, zR, controlparameters, createPlot)
    
    controlparameters.Rstart = radiusS;
    epsilon = 1e-10;
    if thetaR - thetaS > 180
        % Generate true response
        [~, tfmag, ~, fvec, tfcomplex] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaR, radiusR, radiusS, zS, zR, controlparameters, createPlot);
        
        % Generate reference boundary responses
        [~, tfmagDiffRef, ~, ~, ~] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaS + 180 + epsilon, radiusR, radiusS, zS, zR, controlparameters, createPlot);
        [~, tfmagDirRef, ~, ~, ~] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaS + 180 - epsilon, radiusR, radiusS, zS, zR, controlparameters, createPlot);
        
        % Create scaled response
        DiffRef = tfmagDiffRef.diff1;
        DirRef = tfmagDirRef.direct;
        shift = DirRef - DiffRef;
        scaledResponse = tfmag.diff1 + shift;
        truth = 10;

        % Interpolate between the true and scaled responses
        i = min(1, (thetaR - thetaS - 180) / (wedgeIndex - thetaS - 180 - truth));
        tfmag = (1 - i) * scaledResponse + i * tfmag.diff1;
        
        % Adjust magnitude of tfcomplex
        phase = angle(tfcomplex.diff1);
        tfcomplex = PolarToComplex(10.^(tfmag / 20), phase);
    else
        % Remove diffracted sound if not in the shadow zone 
        [~, tfmag, ~, fvec, tfcomplex] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaR, radiusR, radiusS, zS, zR, controlparameters, createPlot);
        tfcomplex.diff1 = zeros(size(tfcomplex.diff1));
        tfmag.diff1 = mag2db(abs(tfcomplex.diff1));
    end
end
