function output = SumAudioOut(audio, fields)

    numFields = length(fields);
    field = fields{1};
    output = [audio.(field).L audio.(field).R];
    for i = 2:numFields
        field = fields{i};
        output = output + [audio.(field).L audio.(field).R];
    end
end