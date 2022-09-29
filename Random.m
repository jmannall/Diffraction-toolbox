function output = Random(range, numOutputs)
    check = range(2) - range(1);
    if check < 0
        error('Number in the first index position must be smaller than the number in the second index position');
    end
    % Create a random set of 'numOutputs' numbers for the given range
    output = range(:,1) +  (range(:,2) - range(:,1)) .* rand(numOutputs,1);
end