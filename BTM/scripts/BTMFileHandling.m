function [inFilePath, fileName, savePath, loadPath, resultExists] = BTMFileHandling(file, index)

    % Create result directory if it doesn't exist
    CheckFileDir('results')

    % Create file info
    [inFilePath,fileName] = fileparts(file);
    fileStem = [fileName, '_', num2str(index)];
    savePath = ['results\', fileName, filesep, fileStem];
    loadPath = ['results\', fileName, filesep, fileStem, '.mat'];

    % Create save directory if doesn't exist
    CheckFileDir(['results\', fileName])

    % Check if previous result exists
    resultExists = isfile(loadPath);
end