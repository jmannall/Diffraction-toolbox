function output = FuncF(X) % 14 or 18 - av 16?

    numFreq = length(X);
    output = zeros(1,numFreq);
    for i = 1:numFreq
        if X(i) < 0.8 % 1
            output(i) = sqrt(pi * X(i)) * (1 - sqrt(X(i)) / (0.7 * sqrt(X(i)) + 1.2)); % 8
        else
            output(i) = 1 - 0.8 / ((X(i) + 1.25)^2); % 4
        end
    end
    output = output .* (exp(1i * (pi / 4) .* (1 - sqrt(X) ./ (X + 1.4)))); % 9
end