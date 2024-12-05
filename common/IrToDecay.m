%% Calculate the decay curve of an impulse response

function [decay, energy] = IrToDecay(ir)  

    structInput = isstruct(ir);

    if structInput
        fields = fieldnames(ir);
        numFields = length(fields);
        numReceivers = size(ir.(fields{1}), 2);
        for n = 1:numFields
            field = fields{n};
            [decay.(field), energy.(field)] = CalcDecay(ir.(field));
        end
    else
        [decay, energy] = CalcDecay(ir);
    end
end

function [E, energy] = CalcDecay(ir)
    temp = cumtrapz(ir(end:-1:1,:).^2); % decay curve
    energy = max(temp);
    z = temp(end:-1:1,:);
    E = 10.*log10(z); % put into dB
    E = E-max(E); % normalise to max 0
    if any(isinf(E))
        E = E(1:find(isinf(E),1,'first')-1,:); % remove trailing infinite values
    end
end