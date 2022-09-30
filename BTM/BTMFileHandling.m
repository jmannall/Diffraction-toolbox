function [inFilePath, fileName, savePath, loadPath, resultExists] = BTMFileHandling(file, index)

    if ~exist('results', 'dir')
           mkdir results
    end

    [inFilePath,fileName] = fileparts(file);
    fileStem = [fileName, '_', num2str(index)];
    savePath = ['results\', fileName, filesep, fileStem];
    loadPath = ['results\', fileName, filesep, fileStem, '.mat'];

    if ~exist(['results\', fileName], 'dir')
           mkdir(['results\' fileName])
    end

    output = exist(loadPath, "file");
    if output == 2
        resultExists = true;
    else
        resultExists = false;
    end
end