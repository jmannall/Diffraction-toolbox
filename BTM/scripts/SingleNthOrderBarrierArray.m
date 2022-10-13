function [result, geometry] = SingleNthOrderBarrierArray(geometry, barrierHeight, thetaS, thetaR, zS, zR, controlparameters)
    
    numInputs = length(geometry.barrierRadius);
        
    % Create file info
    mFile = mfilename('fullpath');
    index = DataHash({geometry, barrierHeight, thetaS, thetaR, zS, zR, controlparameters});
    [fileName, savePath, loadPath, resultExists, filesPerSave, numSaves, extra, rtemplate] = BTMArrayFileHandling(mFile, index, numInputs);

    saveCount = 1;    
    if resultExists
        load(loadPath, "result");
        result = UniformResult(result);
        saveCount = CreateSaveCount(result, numInputs, numSaves, filesPerSave);
    end
    if saveCount <= numSaves
        barrierRadiusi = ReshapeForParfor(geometry.barrierRadius, extra, filesPerSave);
        radiusSi = ReshapeForParfor(geometry.radiusS, extra, filesPerSave);
        radiusRi = ReshapeForParfor(geometry.radiusR, extra, filesPerSave);
        tic
        parfor j = saveCount:numSaves
            result = repmat(rtemplate, 1, 1);
            count = (j - 1) * filesPerSave;
            start = rem(count, filesPerSave);
            batchSize = min(filesPerSave, numInputs - count);
            for i = (start + 1):batchSize
                barrierRadius = barrierRadiusi(i,j);
                radiusS = radiusSi(i,j);
                radiusR = radiusRi(i,j);

                [result(i).ir, result(i).tfmag, result(i).tvec, result(i).fvec, result(i).tfcomplex] = SingleNthOrderBarrier(barrierHeight,barrierRadius,thetaS,thetaR,radiusS,radiusR,zS,zR,controlparameters, false);
                result(i).i = (j - 1) * filesPerSave + i;
            end
            savepathI = [savePath, '_', num2str(j)];
            ParforSaveArray(savepathI, result, geometry);
        end
        toc
        result = ProcessArrayResults(fileName, index, savePath, numSaves, geometry);
    end 
end