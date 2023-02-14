function audio = ConvolveStereoIR(input, ir, windowLength)
    index = DataHash({input, ir, windowLength});
    file = mfilename('fullpath');
    [~,fileName] = fileparts(file);
    fileStem = [fileName, '_', num2str(index)];
    savePath = ['results', filesep, fileName, filesep, fileStem];
    loadPath = ['results', filesep, fileName, filesep, fileStem, '.mat'];
    
    CheckFileDir(['results', filesep, fileName])

    output = exist([cd filesep loadPath], "file");

    if output ~=2
        audio.L = ConvolveIR(input, ir.L, windowLength);
        audio.R = ConvolveIR(input, ir.R, windowLength);
        disp('Save audio')
        save(savePath, "audio")
    else
        disp('Load audio')
        load(loadPath, 'audio')
    end
end