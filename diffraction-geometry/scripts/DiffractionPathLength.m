function pathLength = DiffractionPathLength(receivers, source, zA)

    numReceivers = length(receivers);
    const = zeros(numReceivers, 2);
    apex = [const zA];
    pathLength = PathLength(receivers, apex) + PathLength(apex, source);
end