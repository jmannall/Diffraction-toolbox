function scale = CalculateInterpolationFactor(gParameters)
    theta = [gParameters.thetaS + gParameters.thetaR, gParameters.thetaR - gParameters.thetaS];

    scale = sum((sign(theta - pi) ./ abs(cos(gParameters.v * pi) - cos(gParameters.v * theta)))) ^ 2;
end