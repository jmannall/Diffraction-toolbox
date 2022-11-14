function [azimuth, elevation] = CalculateAzimuthElevation(direction, position, target)
        
    % Calculate horizontal rotation to find azimuth
    z = direction(3);
    positionProj = [position(1:2) z];
    targetProj = [target(1:2) z];
    dirVector = positionProj - targetProj;
    magDirVector = norm(dirVector);
    if magDirVector == 0
        azimuth = 0;
    else
        dirVector = dirVector / norm(dirVector);
        aDot = dot(direction, dirVector);
        aSign = cross(direction, dirVector);
        azimuth = 180 - acosd(aDot);
        
        if aSign(3) > 0
            azimuth = 360 - azimuth;
        end
    end
    
    % Calculate vertical rotation to find elevation
    x = direction(1:2);
    positionProj = [x position(3)];
    targetProj = [x target(3)];
    dirVector = positionProj - targetProj;
    magDirVector = norm(dirVector);
    if magDirVector == 0
        elevation = 0;
    else    
        dirVector = dirVector / norm(dirVector);
        aDot = dot(direction, dirVector);
        aSign = cross(direction, dirVector);
        elevation = 180 - acosd(aDot);
        
        if aSign(1) > 0
            elevation = 360 - elevation;
        end
    end
    
    figure
    plot3(target(:,1), target(:,2), target(:,3), 'o')
    hold on
    plot3(position(:,1), position(:,2), position(:,3), 'o')
    plot3(position(:,1) + [0 direction(1)], position(:,2) + [0 direction(2)], position(:,3) + [0 direction(3)])
    plot3([position(:,1) target(:,1)], [position(:,2) target(:,2)], [position(:,3) target(:,3)])
end