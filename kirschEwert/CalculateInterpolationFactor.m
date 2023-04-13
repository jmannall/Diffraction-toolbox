function scale = CalculateInterpolationFactor(gParameters) % 21
    theta = [gParameters.thetaS + gParameters.thetaR, gParameters.thetaR - gParameters.thetaS];

    scale = sum((sign(theta - pi) ./ abs(cos(gParameters.v * pi) - cos(gParameters.v * theta)))) ^ 2; % 10 * 2 + 1
end