function input = NNInputFromGeometry(wedgeIndex, wedgeLength, thetaR, thetaS, rS, rR, zS, zR, doReflection)
    bendingAngle = thetaR - thetaS;
    minAngle = min(thetaS, wedgeIndex - thetaR);
    
    if length(rS) == 1
        const = ones(size(rR));
        rS = rS * const;
        zS = zS * const;
    end
    % r1 defined as min(rS, rR)
    sourceIsROne = rS < rR;
    [rOne, rTwo, zOne, zTwo] = deal(zeros(size(sourceIsROne)));
    if sum(sourceIsROne) == 0
        rOne = rR;
        rTwo = rS;
        zOne = zR;
        zTwo = zS;
    elseif sum(sourceIsROne) == length(sourceIsROne)
        rOne = rS;
        rTwo = rR;
        zOne = zS;
        zTwo = zR;
    else
        rOne(sourceIsROne, 1) = rS(sourceIsROne);
        rOne(~sourceIsROne, 1) = rR(~sourceIsROne);
        rTwo(sourceIsROne, 1) = rR(sourceIsROne);
        rTwo(~sourceIsROne, 1) = rS(~sourceIsROne);

        % Includes floor reflection (+ wedgeLength)
        zOne(sourceIsROne, 1) = zS(sourceIsROne, 1);
        zOne(~sourceIsROne, 1) = zR(~sourceIsROne, 1);
        zTwo(sourceIsROne, 1) = zR(sourceIsROne, 1);
        zTwo(~sourceIsROne, 1) = zS(~sourceIsROne, 1);
    end
   
    if doReflection
        zOne = zOne + wedgeLength;
        zTwo = zTwo + wedgeLength;
        longWedgeLength = 2 * wedgeLength;
    else
        longWedgeLength = wedgeLength;
    end
    
    % Select other corner if z1 greater than half the wedgeLength (i.e
    % nearer the other corner)
    flipCorner = zOne > longWedgeLength / 2;
    if sum(flipCorner) > 1
        zOne(flipCorner, 1) = longWedgeLength - zOne(flipCorner, 1);
        zTwo(flipCorner, 1) = longWedgeLength - zTwo(flipCorner, 1);
    end
    numReceivers = length(thetaR);
    const = ones(numReceivers, 1);
    input = ([deg2rad(wedgeIndex) .* const, deg2rad(bendingAngle), deg2rad(minAngle), longWedgeLength .* const, rOne, rTwo, zOne, zTwo])';
    input = dlarray(single(input), "CB");
end