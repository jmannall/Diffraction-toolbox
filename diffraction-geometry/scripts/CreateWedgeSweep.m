function [result, geometry, pathLength, validPath] = CreateWedgeSweep(wedgeIndex, minAngle, minBendingAngle, rS, rR, controlparameters, n, interpolated)
    
    bendingAngle = linspace(minBendingAngle, wedgeIndex - minAngle - 0.001, n);
    geometry = GeometryWedge(wedgeIndex, bendingAngle, minAngle, false, false);

    disp('Create wedge sweep')
    wedgeLength = 20;
    [zS, zR] = deal(wedgeLength / 2);
    if interpolated
        [result, geometry] = SingleWedgeInterpolatedArray(geometry, wedgeLength, rS, rR, zS, zR, controlparameters);
    else
        [result, geometry] = SingleWedgeArray(geometry, wedgeLength, rS, rR, zS, zR, controlparameters);
    end

    source = [rS * sind(minAngle); rS * cosd(minAngle)];
    receiver = [rR * sind(minAngle + bendingAngle); rR * cosd(minAngle + bendingAngle)];

    pathLength.dir = vecnorm(source - receiver);
    pathLength.spec = vecnorm([-source(1); source(2)] - receiver);
    pathLength.diff = (rS + rR) * ones(1, n);

    validPath.dir = bendingAngle <= 180;
    validPath.spec = bendingAngle <= 180 - 2 * minAngle;
    validPath.diff = bendingAngle > 0;
    validPath.diffShadow = bendingAngle > 180;
end