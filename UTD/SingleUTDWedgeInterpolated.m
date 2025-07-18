function [tfmag, test, fvec, tfcomplex] = SingleUTDWedgeInterpolated(thetaS, thetaR, radiusS, radiusR, wedgeIndex, phii, controlparameters, A, B)
    
    input = [1; zeros(5, 1)];
    c = controlparameters.c;
    nfft = controlparameters.nfft;
    fs = controlparameters.fs;
    epsilon = 1e-10;
    if isfield(controlparameters, 'fvec')
        fvec = controlparameters.fvec;
    else
        fvec = [125 500 2000 11050];
    end
    k = 2 * pi * fvec/ c; % Precompute
    if thetaR - thetaS > 180
        % Generate true response
        [tfmag, fvec, tfcomplex] = SingleUTDWedge(thetaS, thetaR, radiusS, radiusR, wedgeIndex, phii, controlparameters);
        
        % Generate reference boundary responses
        [DiffRef, ~, ~] = SingleUTDWedge(thetaS, thetaS + 180 + epsilon, radiusS, radiusR, wedgeIndex, phii, controlparameters);

        l = radiusS / sind(phii);
        m = radiusR / sind(phii);
        pathLength = l + m;
        %pathLength = sqrt((radiusS + radiusR) ^ 2 + abs(cosd(phii) * (radiusR - radiusS)) ^ 2);
        %DirRef = mag2db(abs(exp(-1i * k * pathLength) / pathLength));
        DirRef = ones(size(k)) .* mag2db(abs(1 / pathLength));
        %[~, tfmagDirRef, ~, ~, ~] = SingleWedge(wedgeLength, wedgeIndex, thetaS, thetaS + 180 - epsilon, radiusR, radiusS, zS, zR, controlparameters, createPlot);
        
        % Create scaled response
        shift = DirRef - DiffRef;
        scaledResponse = tfmag + shift;

        trueRes = db2mag(tfmag) * pathLength;
        % Interpolate between the true and scaled responses
        i = (thetaR - thetaS - 180) / (wedgeIndex - thetaS - 180);
        tfmag = ((1 - i) * scaledResponse + i * tfmag);
        
        sbRes = db2mag(DiffRef) * pathLength;


        scaleRes = mag2db(trueRes ./ (sbRes * pathLength));

        scaleRes = trueRes ./ sbRes;
        interpResNew = scaleRes .^ (1 - i) .* trueRes .^ i;
        interpRes = (1 - i) * scaleRes + i * trueRes;
        test = mag2db(interpResNew / pathLength);
        % Adjust magnitude of tfcomplex
        phase = angle(tfcomplex);
        tfcomplex = PolarToComplex(10.^(tfmag / 20), phase);
    else
        % Direct sound if not in the shadow zone
        l = radiusS * sind(phii);
        m = radiusR * sind(phii);

        pathLength = sqrt(l ^ 2 + m ^ 2 - 2 * l * m * cosd(thetaR - thetaS));
        tfcomplex = exp(-1i * k * pathLength) / pathLength;
        tfmag = mag2db(abs(tfcomplex));
    end
end
