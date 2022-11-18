function azimuthDeg = CalculateAzimuth(xInput, yInput)
    
    numPositions = length(xInput);
    azimuthDeg = zeros(numPositions, 1);
    for i = 1:numPositions
        x = xInput(i);
        y = yInput(i);

        if x ~= 0
            arg = y ./ x;
            if x < 0
                azimuth = pi + atan(arg);
            else % x > 0
                if y > 0
                    azimuth = atan(arg);
                elseif y < 0
                    azimuth = 2 * pi + atan(arg);
                else
                    azimuth = 0;
                end
            end
        else
            if y > 0
                azimuth = pi / 2;
            elseif y < 0
                azimuth = 3 * pi / 2;
            else % y == 0
                azimuth = 0;
            end
        end
        azimuthDeg(i) = rad2deg(azimuth);
    end
end