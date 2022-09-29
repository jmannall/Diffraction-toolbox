function x = ZPKv2Initialise(i, startingValues)
    x = startingValues;
    if i > 1
        noise = unifrnd(-0.2, 0.2, size(x));
        x = x + noise;
    end
end