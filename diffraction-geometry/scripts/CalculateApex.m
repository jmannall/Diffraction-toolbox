function zA = CalculateApex(radiusS, radiusR, zS, zR, zE)
    if zS < 0 && zR < 0
        zA = 0;
    elseif zS > zE && zR > zE
        zA = zE;
    else
        dZ = abs(zR - zS);
        dZA = dZ * radiusS / (radiusS + radiusR);
        if zS > zR
            dZA = -dZA;
        end
        zA = min(zE, max(0, zS + dZA));
    end
end