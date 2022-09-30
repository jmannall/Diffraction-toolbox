function result = ProcessArrayResults(fileName, index, savePath, numSaves, geometry)
    
    % Compile saves into a single file
    resultAll = [];
    disp('Compiling saves. DO NOT STOP SCRIPT!')
    loadStem = ['results\', fileName, filesep, fileName, '_', num2str(index), '_'];
    for i = 1:numSaves
        loadpathI = [loadStem, num2str(i), '.mat'];
        load(loadpathI, "result");
        resultAll = [resultAll, result];
    end
    result = resultAll;

    % Save results
    save(savePath, "result", "geometry", '-v7.3')
    disp('Result saved')

    % Clear up and delete files
    for i = 1:numSaves
        loadpathI = [loadStem, num2str(i), '.mat'];
        delete(loadpathI);
    end
end