function output = RandomLoguniformDistribution(input, numOutputs)

    check = input(2) - input(1);
    if min(check) <= 0
        error('Number in the first index position must be smaller than the number in the second index position');
    end
    % Create a log uniformly distributed random set of 'numOutputs' numbers for the given range
    lower = input(1);
    upper = input(2);
    % pd = makedist('Loguniform', 'Lower', lower, 'Upper', upper);
    pd = makedist('Uniform', 'Lower', 0, 'Upper', 1);
    y = random(pd, numOutputs, 1);
    output = 10 .^ (y * (log10(upper) - log10(lower)) + log10(lower));
end