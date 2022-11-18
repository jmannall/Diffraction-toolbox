function [zA, phii] = CalculateApex(radiusS, radiusR, zS, zR, zE)
    
    numPositions = length(zR);
    zA = zeros(numPositions, 1);
    for i = 1:numPositions
        if zS < 0 && zR(i) < 0
            zA(i) = 0;
        elseif zS > zE && zR(i) > zE
            zA(i) = zE;
        else
            dZ = abs(zR(i) - zS);
            dZA = dZ * radiusS / (radiusS + radiusR(i));
            if zS > zR(i)
                dZA = -dZA(i);
            end
            zA(i) = min(zE, max(0, zS + dZA));
        end
    end
    dZ = abs(zR - zA);
    phii = atand(radiusR ./ dZ);
end