function output = RandomTriangularDistribution(input, increasing, num)
    
    check = input(:,2) - input(:,1);
    if min(check) <= 0
        error('Number in the first index position must be smaller than the number in the second index position');
    end
    a = 0;
    c = 1;
    if increasing
        b = c;
    else
        b = a;
    end
    pd = makedist('Triangular', 'A', a, 'B', b, 'C', c);
    % Create a triangularly distributed random set of 'numOutputs' numbers for the given range
    output = input(:,1) + (input(:,2) - input(:,1)) .* random(pd, num, 1);
end