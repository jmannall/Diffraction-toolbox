function [geometry, ginput] = NNGeometrySimple(gParameters, numObservations)
    constant = ones(numObservations,1);
    epsilon = gParameters.epsilon;

    % Create wedges
    rangeWedgeLength = gParameters.wedgeLength;
    deltaWedgeLength = rangeWedgeLength(2) - rangeWedgeLength(1);
    rangeWedgeIndex = gParameters.wedgeIndex;
    wedgeLength = Random(rangeWedgeLength, numObservations);
    % Transform wedgeLength to weight lower values
    wedgeLength = -deltaWedgeLength * sqrt(1 - ((wedgeLength - rangeWedgeLength(1)) / deltaWedgeLength).^2) + rangeWedgeLength(2);
    wedgeIndex = Random(rangeWedgeIndex, numObservations);
    
    rangeThetaS = [gParameters.thetaS * constant, wedgeIndex - epsilon];
    rangeThetaR = [gParameters.thetaR * constant, wedgeIndex - epsilon];
    thetaS = Random(rangeThetaS, numObservations);
    thetaR = Random(rangeThetaR, numObservations);
    
    rangeRadiusS = gParameters.radiusS;
    rangeRadiusR = gParameters.radiusR;
    radiusS = Random(rangeRadiusS, numObservations);
    radiusR = Random(rangeRadiusR, numObservations);
    
    rangeZS = [-gParameters.zS * constant, wedgeLength + gParameters.zS];
    rangeZR = [-gParameters.zR * constant, wedgeLength + gParameters.zR];
    %rangeZS = [0 * constant, wedgeLength];
    %rangeZR = [0 * constant, wedgeLength];
    zS = Random(rangeZS, numObservations);
    zR = Random(rangeZR, numObservations);

    minAngle = min(thetaS, thetaR);
    bendingAngle = max(thetaS, thetaR) - minAngle;
    midPoint = wedgeLength / 2;
    z1 = abs(midPoint - zS);
    z2 = abs(midPoint - zR);
    deltaZ = abs(zR - zS);

    check = bendingAngle == 180;
    if sum(check) > 0
        thetaR(check) = thetaR(check) + epsilon;
        minAngle = min(thetaS, thetaR);
        bendingAngle = max(thetaS, thetaR) - minAngle;
    end
    check = thetaS == 180;
    if sum(check) > 0
        thetaS(check) = thetaS(check) + epsilon;
        minAngle = min(thetaS, thetaR);
        bendingAngle = max(thetaS, thetaR) - minAngle;
    end
    check = wedgeIndex - thetaS == 180;
    if sum(check) > 0
        thetaS(check) = thetaS(check) + epsilon;
        minAngle = min(thetaS, thetaR);
        bendingAngle = max(thetaS, thetaR) - minAngle;
    end
    check = thetaR == 180;
    if sum(check) > 0
        thetaR(check) = thetaR(check) + epsilon;
        minAngle = min(thetaS, thetaR);
        bendingAngle = max(thetaS, thetaR) - minAngle;
    end
    check = wedgeIndex - thetaR == 180;
    if sum(check) > 0
        thetaR(check) = thetaR(check) + epsilon;
        minAngle = min(thetaS, thetaR);
        bendingAngle = max(thetaS, thetaR) - minAngle;
    end

    geometry = struct('wedgeLength', wedgeLength, 'wedgeIndex', wedgeIndex, ...
        'thetaS', thetaS, 'thetaR', thetaR, 'radiusS', radiusS, 'radiusR', radiusR, 'zS', zS, 'zR', zR);
    
    ginput = struct('wedgeLength', wedgeLength, 'wedgeIndex', wedgeIndex, ...
        'minAngle', minAngle, 'bendingAngle', bendingAngle, 'z1', z1, 'z2', z2, 'deltaZ', deltaZ);
end