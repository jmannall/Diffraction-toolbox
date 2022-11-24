function [fileName, savePath, loadPath, resultExists, filesPerSave, numSaves, extra, rtemplate] = BTMArrayFileHandling(file, index, numInputs)

    % Create result directory if it doesn't exist
    CheckFileDir('results')

    % Create file info
    [~,fileName] = fileparts(file);
    fileStem = [fileName, '_', num2str(index)];
    savePath = ['results', filesep, fileName, filesep, fileStem];
    loadPath = ['results', filesep, fileName, filesep, fileStem, '.mat'];

    % Create save directory if doesn't exist
    CheckFileDir(['results', filesep, fileName])

    % Check if previous complete result exists
    output = exist([cd filesep loadPath], "file");

    % Loop through to check for previous unifinished attempts
    i = 0;
    extra = 0;
    filesPerSave = numInputs / 20;
    numSaves = ceil(numInputs / filesPerSave);
    if output ~= 2
        test = true;
        loadPath = loadPath(1:end-4);
        while test == 2 && i <= numSaves
            i = i + 1;
            loadpathI = [loadPath, '_', num2str(i), '.mat'];
            test = exist([cd filesep loadpathI], "file");
        end
        % Check if previous result exists
        loadPath = [loadPath, '_', num2str(max(1,i - 1)), '.mat'];
        output = exist([cd filesep loadPath], "file");
        extra = filesPerSave - rem(numInputs, filesPerSave);
    end

    if output == 2
        resultExists = true;
    else
        resultExists = false;
    end

    % Create result template
    rtemplate.ir = [];
    rtemplate.tfmag = [];
    rtemplate.tvec = [];
    rtemplate.fvec = [];
    rtemplate.tfcomplex = [];
    rtemplate.i = [];
end