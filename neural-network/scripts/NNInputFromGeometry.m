function input = NNInputFromGeometry(wedgeIndex, wedgeLength, thetaR, thetaS, rS, rR, zS, zR)
    bendingAngle = thetaR - thetaS;
    minAngle = min(thetaS, wedgeIndex - thetaR);
    
    % r1 defined as min(rS, rR)
    sourceIsROne = rS < rR;
    rOne(sourceIsROne, 1) = rS;
    rOne(~sourceIsROne, 1) = rR(~sourceIsROne);
    rTwo(sourceIsROne, 1) = rR(sourceIsROne);
    rTwo(~sourceIsROne, 1) = rS;
    
    % Includes floor reflection (+ wedgeLength)
    zOne(sourceIsROne, 1) = zS + wedgeLength;
    zOne(~sourceIsROne, 1) = zR(~sourceIsROne, 1) + wedgeLength;
    zTwo(sourceIsROne, 1) = zR(sourceIsROne, 1) + wedgeLength;
    zTwo(~sourceIsROne, 1) = zS + wedgeLength;
    
    longWedgeLength = 2 * wedgeLength;
    
    % Select other corner if z1 greater than half the wedgeLength (i.e
    % neraer the other corner)
    flipCorner = zOne > wedgeLength;
    zOne(flipCorner, 1) = longWedgeLength - zOne(flipCorner, 1);
    zTwo(flipCorner, 1) = longWedgeLength - zTwo(flipCorner, 1);
    
    numReceivers = length(thetaR);
    const = ones(numReceivers, 1);
    input = ([deg2rad(wedgeIndex) * const, deg2rad(bendingAngle), deg2rad(minAngle), longWedgeLength * const, rOne, rTwo, zOne, zTwo])';
    input = dlarray(single(input), "CB");
end