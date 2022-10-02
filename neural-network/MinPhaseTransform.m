%% Constrain a zero or pole to lie within the unit circle. Traceable by dlgradient

function [realOutput, imagOutput] = MinPhaseTransform(realPart, imagPart)
    epsilon = 1e-8;
    mag = sqrt(realPart.^2 + imagPart.^2);
    magEps = sqrt((realPart + epsilon).^2 + imagPart.^2);

    realOutput = (1 - epsilon) .* realPart .* tanh(mag) ./ magEps;
    imagOutput = (1 - epsilon) .* imagPart .* tanh(mag) ./ magEps;
end