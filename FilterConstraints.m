function x = FilterConstraints(x)
    
    % Set filter parameter order as zero, pole, zero, pole
    options = sort(x(1:4));
    x(1:4) = options;

    % Set k to be positive
    x(5) = abs(x(5));

    % Apply Lower and Upper Bound Limits
    x = max(x, -1);
    x = min(x, 1);
end