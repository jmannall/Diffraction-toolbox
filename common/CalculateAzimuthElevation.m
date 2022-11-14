function [azimuth, elevation] = CalculateAzimuthElevation(direction, position, target)
    numPositions = size(direction, 1);
    const = ones(numPositions, 2);

    % Calculate horizontal rotation to find azimuth
    z = direction(:,3);
    positionProj = [position(:,1:2) z];
    targetProj = [const .* target(1:2) z];
    dirVector = positionProj - targetProj;
    magDirVector = vecnorm(dirVector, 2, 2);
    dirVector = dirVector ./ magDirVector;
    aDot = dot(direction, dirVector, 2);
    aSign = cross(direction, dirVector, 2);
    azimuth = 180 - acosd(aDot);
    
    azimuth(aSign(:,3) > 0) = 360 - azimuth(aSign(:,3) > 0);
    azimuth(magDirVector == 0) = 0;
    
    const = ones(numPositions, 1);

    % Calculate vertical rotation to find elevation
    x = direction(:, 1:2);
    positionProj = [x position(:,3)];
    targetProj = [x const .* target(3)];
    dirVector = positionProj - targetProj;
    magDirVector = vecnorm(dirVector, 2, 2);  
    dirVector = dirVector ./ magDirVector;
    aDot = dot(direction, dirVector, 2);
    aSign = cross(direction, dirVector, 2);
    elevation = 180 - acosd(aDot);

    elevation(aSign(:,1) > 0) = 360 - elevation(aSign(:,1) > 0);
    elevation(magDirVector == 0) = 0;
end