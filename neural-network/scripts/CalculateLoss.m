function loss = CalculateLoss(prediction, target, filter)

    if nargin < 3
        filter = true(1, size(target, 2));
    end
    structInput = isstruct(prediction);
    if structInput
        fields = fieldnames(prediction);
        numFields = length(fields);
        for i = 1:numFields
            idx = fields{i};
            loss.if.(idx) = abs(target(:,filter) - prediction.(idx)(:,filter));
            loss.i.(idx) = mean(loss.if.(idx), 1);
            loss.f.(idx) = mean(loss.if.(idx), 2);
            loss.mean.(idx) = mean(loss.if.(idx), 'all');
            loss.ordered.(idx) = sort(loss.i.(idx));
        end
    else
        loss.if = abs(target(:,filter) - prediction(:,filter));
        loss.i = mean(loss.if, 1);
        loss.f = mean(loss.if, 2);
        loss.mean = mean(loss.if, 'all');
    end
end