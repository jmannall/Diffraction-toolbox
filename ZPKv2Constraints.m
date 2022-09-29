function x = ZPKv2Constraints(x, input)
    x(1:4) = max(-0.99, min(x(1:4), 0.99));
    x(5) = max(0, min(x(5), 1));
    store = sort(x(1:4));
    x(1:4) = fliplr(store);
    prediction = x * input;
    checkMax = max(prediction, [], 2);
    idxMax = find(checkMax > 0.99);
    checkMin = min(prediction, [], 2);
    idxMin = find(checkMax < -0.99);
    if sum(idxMax) > 0
        scale = 1 ./ checkMax(idxMax);
        x(idxMax,:) = scale .* x(idxMax,:);
    end
    if sum(idxMin) > 0
        scale = -1 ./ checkMin(idxMin);
        x(idxMin,:) = scale .* x(idxMin,:);
    end

    prediction = x * input;
    checkMax = max(prediction, [], 2);
    checkMin = min(prediction, [], 2);
end