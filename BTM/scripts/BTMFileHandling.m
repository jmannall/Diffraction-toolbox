function [inFilePath, fileName, savePath, loadPath, resultExists] = BTMFileHandling(file, index)

    % Create result directory if it doesn't exist
    if ~exist('results', 'dir')
           mkdir results
    end

    % Create file info
    [inFilePath,fileName] = fileparts(file);
    fileStem = [fileName, '_', num2str(index)];
    savePath = ['results\', fileName, filesep, fileStem];
    loadPath = ['results\', fileName, filesep, fileStem, '.mat'];

    % Create save directory if doesn't exist
    if ~exist(['results\', fileName], 'dir')
           mkdir(['results\' fileName])
    end

    % Check if previous result exists
    output = exist(loadPath, "file");
    if output == 2
        resultExists = true;
    else
        resultExists = false;
    end
end