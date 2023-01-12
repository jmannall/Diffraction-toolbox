function [result, geometry] = SingleWedgeInterpolatedArray(geometry, wedgeLength, radiusS, radiusR, zS, zR, controlparameters)
    
    numInputs = length(geometry.wedgeIndex);
            
    % Create file info
    mFile = mfilename('fullpath');
    index = DataHash({geometry, wedgeLength, radiusS, radiusR, zS, zR, controlparameters});
    [fileName, savePath, loadPath, resultExists, filesPerSave, numSaves, extra, rtemplate] = BTMArrayFileHandling(mFile, index, numInputs);

    saveCount = 1;
    if resultExists
        load(loadPath, "result");
        %result = UniformResult(result);
        saveCount = CreateSaveCount(result, numInputs, numSaves, filesPerSave);
    end
    if saveCount <= numSaves
        wedgeIndexi = ReshapeForParfor(geometry.wedgeIndex, extra, filesPerSave);
        thetaSi = ReshapeForParfor(geometry.source, extra, filesPerSave);
        thetaRi = ReshapeForParfor(geometry.receiver, extra, filesPerSave);
        tic
        for j = saveCount:numSaves
            result = repmat(rtemplate, 1, 1);
            count = (j - 1) * filesPerSave;
            start = rem(count, filesPerSave);
            batchSize = min(filesPerSave, numInputs - count);
            for i = (start + 1):batchSize
                wedgeIndex = wedgeIndexi(i,j);
                thetaS = thetaSi(i,j);
                thetaR = thetaRi(i,j);
                [result(i).tfmag, result(i).fvec, result(i).tfcomplex] = SingleWedgeInterpolated(wedgeLength,wedgeIndex,thetaS,thetaR,radiusS,radiusR,zS,zR,controlparameters,false);
                result(i).i = (j - 1) * filesPerSave + i;
            end
            savepathI = [savePath, '_', num2str(j)];
            ParforSaveArray(savepathI, result, geometry);
        end
        toc
        result = ProcessArrayResults(fileName, index, savePath, numSaves, controlparameters, geometry);
    end 
end