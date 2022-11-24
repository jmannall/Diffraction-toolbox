function result = ProcessArrayResults(fileName, index, savePath, numSaves, controlparameters, geometry, NNinput)
    
    % Compile saves into a single file
    resultAll = [];
    disp('Compiling saves. DO NOT STOP SCRIPT!')
    loadStem = ['results',filesep, fileName, filesep, fileName, '_', num2str(index), '_'];
    for i = 1:numSaves
        loadpathI = [loadStem, num2str(i), '.mat'];
        load(loadpathI, "result");
        resultAll = [resultAll, result];
    end
    result = resultAll;

    % Save results
    if controlparameters.saveFiles >= 2
        binary = fliplr(dec2bin(controlparameters.saveFiles));
        if binary(2) == '1'
            if nargin > 6
                save(savePath, "result", "geometry", "NNinput", '-v7.3')
            else
                save(savePath, "result", "geometry", '-v7.3')
            end
            disp('Result saved')
        end
    end
    
    if isfield(result, 'ir')
        if ~isempty([result.ir])
            result = UniformResult(result);
        end
    end
    % Clear up and delete files
    for i = 1:numSaves
        loadpathI = [loadStem, num2str(i), '.mat'];
        delete(loadpathI);
    end
end