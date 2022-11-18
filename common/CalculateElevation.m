function elevationDeg = CalculateElevation(xInput, yInput, zInput)

    numPositions = length(xInput);
    elevationDeg = zeros(numPositions, 1);
    for i = 1:numPositions
        x = xInput(i);
        y = yInput(i);
        z = zInput(i);
        if x == 0 && y == 0
            if z > 0
                elevation = pi / 2;
            elseif z < 0
                elevation = -pi / 2;
            else
                elevation = 0;
            end
        else
            elevation = atan(z ./ sqrt(x.^2 + y.^2));
        end
        elevationDeg(i) = min(90, max(-30, rad2deg(elevation)));
    end
end
