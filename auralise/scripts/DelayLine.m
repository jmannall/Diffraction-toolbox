%% Implement delay line

function [output, ir] = DelayLine(audio, pathLength, windowLength, validPath, c, fs)

    structInput = isstruct(pathLength);

    if structInput
        fields = fieldnames(pathLength);
        numFields = length(fields);
        for n = 1:numFields
            field = fields{n};

            [output.(field), ir.(field)] = RunDelayLine(audio, pathLength.(field), windowLength, validPath.(field), c, fs);
        end
    else
        [output, ir] = RunDelayLine(audio, pathLength, windowLength, validPath, c, fs);
    end
end