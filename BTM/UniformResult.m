%% Make result the same diffraction order and ir lengths across positions

function result = UniformResult(result)

    numResults = length(result);
    [irLength, diffOrder] = deal(zeros(numResults, 1));
    for i = 1:numResults
        diffOrder(i) = length(fieldnames(result(i).ir)) - 3;
        irLength(i) = length(result(i).ir.complete);
    end
    maxIrLength = max(irLength);
    maxDiffOrder = max(diffOrder);
    for i = 1:numResults
        % Add structs for empty diffraction orders
        if diffOrder(i) < maxDiffOrder
            idx = ['diff', num2str(diffOrder(i))];
            paddingIr = zeros(size(result(i).ir.(idx)));
            paddingTfmag = zeros(size(result(i).tfmag.(idx)));
            paddingTfcomplex = zeros(size(result(i).tfcomplex.(idx)));
            for j = diffOrder(i) + 1:maxDiffOrder
                idx = ['diff', num2str(j)];
                result(i).ir.(idx) = paddingIr;
                result(i).tfmag.(idx) = paddingTfmag;
                result(i).tfcomplex.(idx) = paddingTfcomplex;
            end
        end

        % Make all impulse responses the same length
        padding = zeros(maxIrLength - irLength(i), 1);
        result(i).ir.complete = [result(i).ir.complete; padding];
        result(i).ir.direct = [result(i).ir.direct; padding];
        result(i).ir.geom = [result(i).ir.geom; padding];
        for j = 1:maxDiffOrder
            idx = ['diff', num2str(j)];
            result(i).ir.(idx) = [result(i).ir.(idx); padding];
        end
    end
end