function [zA, phii] = CalculateApex(radiusS, radiusR, zS, zR, zE, noClip)
    
    numPositions = length(zR);
    zA = zeros(numPositions, 3);
    if length(zS) == 1
        const = ones(size(zR));
        zS = zS * const;
        zE = zE * const;
        radiusS = radiusS * const;
    end
    if length(zE) == 1
        zE = zE * ones(size(zR));
    end
    if noClip
        for i = 1:numPositions
            dZ = abs(zR(i) - zS(i));
            dZA = dZ * radiusS(i) / (radiusS(i) + radiusR(i));
            if zS(i) > zR(i)
                dZA = -dZA;
            end
            zA(i,3) = zS(i) + dZA;
        end
    else
        for i = 1:numPositions
            if zS(i) < 0 && zR(i) < 0
                zA(i,3) = 0;
            elseif zS(i) > zE(i) && zR(i) > zE(i)
                zA(i,3) = zE(i);
            else
                dZ = abs(zR(i) - zS(i));
                dZA = dZ * radiusS(i) / (radiusS(i) + radiusR(i));
                if zS(i) > zR(i)
                    dZA = -dZA;
                end
                zA(i,3) = min(zE(i), max(0, zS(i) + dZA));
            end
        end
    end
    dZ = abs(zR - zA(:,3)');
    phii = atand(radiusR ./ dZ);
end