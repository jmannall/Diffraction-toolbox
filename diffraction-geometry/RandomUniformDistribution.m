function output = RandomUniformDistribution(input, numOutputs)

    check = input(:,2) - input(:,1);
    if min(check) <= 0
        error('Number in the first index position must be smaller than the number in the second index position');
    end
    % Create a uniformly distributed random set of 'numOutputs' numbers for the given range
    output = input(:,1) +  (input(:,2) - input(:,1)) .* rand(numOutputs,1);
end