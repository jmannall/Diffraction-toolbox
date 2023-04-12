function percentile = CalculatePercentiles(loss, percentiles)

    structInput = isstruct(loss);
    if structInput
        fields = fieldnames(loss);
        numFields = length(fields);
        for i = 1:numFields
            idx = fields{i};
            percentile.(idx) = prctile(loss.(idx), percentiles);
        end
    else
        percentile.(idx) = prctile(loss, percentiles);
    end
end