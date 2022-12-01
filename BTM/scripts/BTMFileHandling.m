function [inFilePath, fileName, savePath, loadPath, resultExists] = BTMFileHandling(file, index)

    % Create result directory if it doesn't exist
    CheckFileDir('results')

    % Create file info
    [inFilePath,fileName] = fileparts(file);
    fileStem = [fileName, '_', num2str(index)];
    savePath = ['results', filesep, fileName, filesep, fileStem];
    loadPath = ['results', filesep, fileName, filesep, fileStem, '.mat'];

    % Create save directory if doesn't exist
    % CheckFileDir(['results', filesep, fileName])

    % Check if previous result exists
    resultExists = isfile(loadPath);
end