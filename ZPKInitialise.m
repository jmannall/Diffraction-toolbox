function x = ZPKInitialise(i, startingValues)
    x = startingValues;
    if i > 1
        noise = unifrnd(-0.001, 0.001, size(x));
        x = x + noise;
    end
end