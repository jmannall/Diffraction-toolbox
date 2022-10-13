function [fileName, savePath, loadPath, resultExists, filesPerSave, numSaves, extra, rtemplate] = BTMArrayFileHandling(file, index, numInputs)

    % Create result directory if it doesn't exist
    if ~exist('results', 'dir')
           mkdir results
    end

    % Create file info
    [~,fileName] = fileparts(file);
    fileStem = [fileName, '_', num2str(index)];
    savePath = ['results\', fileName, filesep, fileStem];
    loadPath = ['results\', fileName, filesep, fileStem, '.mat'];

    % Create save directory if doesn't exist
    if ~exist(['results\', fileName], 'dir')
           mkdir(['results\' fileName])
    end

    % Check if previous complete result exists
    output = exist(loadPath, "file");

    % Loop through to check for previous unifinished attempts
    i = 0;
    extra = 0;
    filesPerSave = 50;
    numSaves = ceil(numInputs / filesPerSave);
    if output ~= 2
        test = 2;
        loadPath = loadPath(1:end-4);
        while test == 2 && i <= numSaves
            i = i + 1;
            loadpathI = [loadPath, '_', num2str(i), '.mat'];
            test = exist(loadpathI, "file");
        end
        loadPath = [loadPath, '_', num2str(max(1,i - 1)), '.mat'];
        output = exist(loadPath, "file");
        extra = filesPerSave - rem(numInputs, filesPerSave);
    end

    % Check if previous result exists
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